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

  size_t num_local_cells = grid->local_cell_glob_indices.size();
  size_t tally_size = num_grps*num_local_cells;
  phi_tally_contrib.resize(tally_size,0.0);
  phi_tally.resize(tally_size,0.0);
  phi_tally_sqr.resize(tally_size,0.0);

  phi_global_initial_value.resize(tally_size,0.0);
  phi_global.resize(tally_size,0.0);
  phi_global_tally_sqr.resize(tally_size,0.0);

  phi_local_relsigma.resize(tally_size,0.0);

  //=================================== Initialize Finite Volume discretization
  chi_log.Log(LOG_0) << "Adding finite volume views.";
  fv_discretization = new SpatialDiscretization_FV;

  fv_discretization->AddViewOfLocalContinuum(grid);

  MPI_Barrier(MPI_COMM_WORLD);
  if (make_pwld)
  {
    chi_log.Log(LOG_0) << "Adding PWL finite element views.";
    //=================================== Initialize pwl discretization
    pwl_discretization = new SpatialDiscretization_PWL;

    pwl_discretization->
      AddViewOfLocalContinuum(grid);

    //=================================== Generate moment wise addresses
    num_moms = 1;
    local_cell_pwl_dof_array_address.resize(num_local_cells,0);
    int block_MG_counter = 0;
    for (size_t lc=0; lc<num_local_cells; lc++)
    {
      local_cell_pwl_dof_array_address[lc] = block_MG_counter;
      auto cell_pwl_view = pwl_discretization->MapFeViewL(lc);

      block_MG_counter += cell_pwl_view->dofs*num_grps*num_moms;
    }

    //=================================== Initialize PWLD tallies
    tally_size = block_MG_counter;
    phi_pwl_tally_contrib.resize(tally_size,0.0);
    phi_pwl_tally.resize(tally_size,0.0);
    phi_pwl_tally_sqr.resize(tally_size,0.0);

    phi_pwl_global.resize(tally_size,0.0);
    phi_pwl_global_tally_sqr.resize(tally_size,0.0);

    phi_pwl_local_relsigma.resize(tally_size,0.0);
  }
  MPI_Barrier(MPI_COMM_WORLD);
}