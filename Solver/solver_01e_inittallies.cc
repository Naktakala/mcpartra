#include "solver_montecarlon.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

#include <unistd.h>

//###################################################################
/**Initializes tallies.*/
void mcpartra::Solver::InitTallies()
{
  chi_log.Log() << "MCParTra: Initializing tallies";

  //=================================== Unknown Manager
  for (int m=0; m<num_moms; ++m)
  {
    uk_man_fv .AddUnknown(chi_math::UnknownType::VECTOR_N, num_grps);
    uk_man_pwld.AddUnknown(chi_math::UnknownType::VECTOR_N, num_grps);
    for (int g=0; g<num_grps; ++g)
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

  //=================================== Initialize tally blocks
  grid_tally_blocks.clear();
  grid_tally_blocks.emplace_back(); //DEFAULT_FVTALLY
  grid_tally_blocks.emplace_back(); //DEFAULT_PWLTALLY
  grid_tally_blocks.emplace_back(); //UNCOLLIDED_FVTALLY
  grid_tally_blocks.emplace_back(); //UNCOLLIDED_PWLTALLY


  //=================================== Initialize Finite Volume discretization
  chi_log.Log(LOG_0) << "Adding finite volume views.";
  fv = SpatialDiscretization_FV::New(grid);

  usleep(1000000);

  //=================================== Tally sizes
  auto fv_tally_size = fv->GetNumLocalDOFs(uk_man_fv);

  grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].Resize(fv_tally_size);
  grid_tally_blocks[TallyMaskIndex[UNCOLLIDED_FVTALLY]].Resize(fv_tally_size);

  MPI_Barrier(MPI_COMM_WORLD);


  //=================================== Initialize pwl discretization
  chi_log.Log(LOG_0) << "Adding PWL finite element views.";
//  pwl = SpatialDiscretization_PWLD::New(grid,
//          chi_math::finite_element::COMPUTE_CELL_MAPPINGS |
//          chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);
  pwl = SpatialDiscretization_PWLD::New(grid,
                                        chi_math::finite_element::COMPUTE_CELL_MAPPINGS);

  usleep(1000000);

  //=================================== Initialize PWLD tallies
  auto fem_tally_size = pwl->GetNumLocalDOFs(uk_man_pwld);
  auto fem_tally_size_global = pwl->GetNumGlobalDOFs(uk_man_pwld);

  chi_log.Log(LOG_0) << "PWL #global-dofs: " << fem_tally_size_global;

  grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].Resize(fem_tally_size);
//  grid_tally_blocks[TallyMaskIndex[UNCOLLIDED_PWLTALLY]].Resize(fem_tally_size);

  MPI_Barrier(MPI_COMM_WORLD);
  chi_log.Log(LOG_0) << "Done sizing tallies";
  usleep(1000000);

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

  chi_log.Log(LOG_0) << "Done initializing tallies.";
  MPI_Barrier(MPI_COMM_WORLD);
}