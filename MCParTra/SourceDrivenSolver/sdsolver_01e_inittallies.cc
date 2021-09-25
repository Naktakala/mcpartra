#include "sdsolver.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include <iomanip>

//###################################################################
/**Initializes tallies.*/
void mcpartra::SourceDrivenSolver::InitTallies()
{
  chi_log.Log() << "MCParTra: Initializing tallies";

  uk_man_fv.unknowns.clear();
  uk_man_pwld.unknowns.clear();
  //=================================== Unknown Managers moments and groups
  for (size_t m=0; m < num_moments; ++m)
  {
    uk_man_fv .AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);
    uk_man_pwld.AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);
    for (size_t g=0; g < num_groups; ++g)
    {
      auto fv_comp_name = std::string("PhiFV");
      auto fem_comp_name = std::string("PhiFEM");

      fv_comp_name += std::string("_m")+std::to_string(m);
      fv_comp_name += std::string("_g")+std::to_string(g);

      fem_comp_name += std::string("_m")+std::to_string(m);
      fem_comp_name += std::string("_g")+std::to_string(g);

      uk_man_fv .SetUnknownComponentTextName(m, g, fv_comp_name);
      uk_man_pwld.SetUnknownComponentTextName(m, g, fem_comp_name);
    }//for g
  }//for m

  //=================================== Unknown Manager groups only
  uk_man_fv_importance.unknowns.clear();
  uk_man_fv_importance.AddUnknown(chi_math::UnknownType::VECTOR_N, num_groups);
  for (size_t g=0; g < num_groups; ++g)
  {
    auto imp_comp_name = std::string("CellImportance");

    imp_comp_name += std::string("_g") + std::to_string(g);

    uk_man_fv_importance.SetUnknownComponentTextName(0, g, imp_comp_name);
  }//for g

  //=================================== Estimating total size requirement
  //                                    for tallies only
  size_t num_local_cells = grid->local_cells.size();
  size_t num_globl_cells = grid->GetGlobalNumberOfCells();

  size_t num_local_nodes = 0;
  for (const auto& cell : grid->local_cells)
    num_local_nodes += cell.vertex_ids.size();

  size_t max_num_local_cells = 0;
  MPI_Allreduce(&num_local_cells,         //sendbuf
                &max_num_local_cells,     //recvbuf
                1,MPI_UNSIGNED_LONG_LONG, //count,datatype
                MPI_MAX,                  //operation
                MPI_COMM_WORLD);          //communicator

  size_t num_globl_nodes = 0;
  MPI_Allreduce(&num_local_nodes,         //sendbuf
                &num_globl_nodes,         //recvbuf
                1,MPI_UNSIGNED_LONG_LONG, //count,datatype
                MPI_SUM,                  //operation
                MPI_COMM_WORLD);          //communicator

  size_t max_num_local_nodes = 0;
  MPI_Allreduce(&num_local_nodes,         //sendbuf
                &max_num_local_nodes,     //recvbuf
                1,MPI_UNSIGNED_LONG_LONG, //count,datatype
                MPI_MAX,                  //operation
                MPI_COMM_WORLD);          //communicator

  double fv_tally_peak_size_estimate =
    static_cast<double>(max_num_local_cells) *
    uk_man_fv.GetTotalUnknownStructureSize() *
    8 *      //bytes per float
    6 *      //vectors per tally
    1.0e-6;  //bytes per megabyte

  double fv_tally_total_size_estimate =
    static_cast<double>(num_globl_cells) *
    uk_man_fv.GetTotalUnknownStructureSize() *
    8 *      //bytes per float
    6 *      //vectors per tally
    1e-6;    //bytes per megabyte

  chi_log.Log() << "MCParTra: Estimated size required by FV tallies: peak "
                << std::fixed << std::setprecision(1)
                << fv_tally_peak_size_estimate << "Mb (total "
                << std::fixed << std::setprecision(1)
                << fv_tally_total_size_estimate << "Mb)";

  if (options.make_pwld)
  {
    double pwl_tally_peak_size_estimate =
      static_cast<double>(max_num_local_nodes) *
      uk_man_pwld.GetTotalUnknownStructureSize() *
      8 *      //bytes per float
      6 /      //vectors per tally
      1000000; //bytes per megabyte

    double pwl_tally_total_size_estimate =
      static_cast<double>(num_globl_nodes) *
      uk_man_pwld.GetTotalUnknownStructureSize() *
      8 *      //bytes per float
      6 /      //vectors per tally
      1000000; //bytes per megabyte

    chi_log.Log() << "MCParTra: Estimated size required by PWL tallies: peak "
                  << std::fixed << std::setprecision(1)
                  << pwl_tally_peak_size_estimate << "Mb (total "
                  << std::fixed << std::setprecision(1)
                  << pwl_tally_total_size_estimate << "Mb)";
  }

  //=================================== Initialize tally blocks
  grid_tally_blocks.clear();
  grid_tally_blocks.emplace_back(); //DEFAULT_FVTALLY
  grid_tally_blocks.emplace_back(); //DEFAULT_PWLTALLY
  grid_tally_blocks.emplace_back(); //UNCOLLIDED_FVTALLY
  grid_tally_blocks.emplace_back(); //UNCOLLIDED_PWLTALLY


  //=================================== Initialize Finite Volume discretization
  chi_log.Log(LOG_0) << "Adding finite volume views.";
  fv = SpatialDiscretization_FV::New(grid);

  //=================================== Tally sizes
  auto fv_tally_size = fv->GetNumLocalDOFs(uk_man_fv);

  grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].Resize(fv_tally_size);

  MPI_Barrier(MPI_COMM_WORLD);


  //=================================== Initialize pwl discretization
  pwl = SpatialDiscretization_PWLD::New(grid,
              chi_math::finite_element::COMPUTE_CELL_MAPPINGS);

  //=================================== Initialize PWLD tallies
  auto fem_tally_size = pwl->GetNumLocalDOFs(uk_man_pwld);
  auto fem_tally_size_global = pwl->GetNumGlobalDOFs(uk_man_pwld);

  grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].Resize(fem_tally_size);

  //=================================== Initialize custom tallies
  auto fv_dof_struct_size = uk_man_fv.GetTotalUnknownStructureSize();
  for (auto& custom_tally : custom_tallies)
  {
    double local_volume = 0.0;
    for (auto& cell : grid->local_cells)
      if (custom_tally.local_cell_tally_mask[cell.local_id])
        local_volume += fv->MapFeView(cell.local_id)->volume;

    double global_volume = 0.0;
    MPI_Allreduce(&local_volume,  //sendbuf
                  &global_volume, //recvbuf
                  1,              //count
                  MPI_DOUBLE,     //datatype
                  MPI_SUM,        //operation
                  MPI_COMM_WORLD);//communicator

    custom_tally.Initialize(fv_dof_struct_size,global_volume);
  }

  chi_log.Log(LOG_0) << "MCParTra: Done initializing tallies.";
}