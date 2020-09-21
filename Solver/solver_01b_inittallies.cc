#include "solver_montecarlon.h"

#include <ChiMesh/MeshHandler/chi_meshhandler.h>

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Initializes tallies.*/
void chi_montecarlon::Solver::InitTallies()
{
  auto handler = chi_mesh::GetCurrentHandler();
  mesh_is_global = handler->volume_mesher->options.mesh_global;

  //=================================== Unknown Manager
  for (int m=0; m<num_moms; ++m)
  {
    uk_man_fv .AddUnknown(chi_math::UnknownType::VECTOR_N, num_grps);
    uk_man_fem.AddUnknown(chi_math::UnknownType::VECTOR_N, num_grps);
    for (int g=0; g<num_grps; ++g)
    {
      auto fv_comp_name = std::string("PhiFV");
      auto fem_comp_name = std::string("PhiFEM");

      fv_comp_name += std::string("_m")+std::to_string(m);
      fv_comp_name += std::string("_g")+std::to_string(g);

      fem_comp_name += std::string("_m")+std::to_string(m);
      fem_comp_name += std::string("_g")+std::to_string(g);

      uk_man_fv .SetUnknownComponentTextName(m, g, fv_comp_name);
      uk_man_fem.SetUnknownComponentTextName(m, g, fem_comp_name);
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
  fv = new SpatialDiscretization_FV;
  fv->AddViewOfLocalContinuum(grid);

  //=================================== Tally sizes
  auto fv_tally_size = fv->GetNumLocalDOFs(grid, &uk_man_fv);

  grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].Resize(fv_tally_size);
  grid_tally_blocks[TallyMaskIndex[UNCOLLIDED_FVTALLY]].Resize(fv_tally_size);

  MPI_Barrier(MPI_COMM_WORLD);

  //=================================== Initialize pwl discretization
  chi_log.Log(LOG_0) << "Adding PWL finite element views.";
  pwl = new SpatialDiscretization_PWL;
  pwl->AddViewOfLocalContinuum(grid);
  pwl->OrderNodesDFEM(grid);

  //=================================== Initialize PWLD tallies
  auto fem_tally_size = pwl->GetNumLocalDOFs(grid, &uk_man_fem);

  chi_log.Log(LOG_0) << "PWL #local-dofs: " << fem_tally_size;

  grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].Resize(fem_tally_size);
  grid_tally_blocks[TallyMaskIndex[UNCOLLIDED_PWLTALLY]].Resize(fem_tally_size);

  chi_log.Log(LOG_0) << "Done initializing tallies.";
  MPI_Barrier(MPI_COMM_WORLD);
}