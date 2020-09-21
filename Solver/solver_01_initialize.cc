#include"solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Initialize the solver.*/
bool chi_montecarlon::Solver::Initialize()
{
  chi_log.Log(LOG_0) << "Initializing MonteCarlo solver.";

  //=================================== Obtain grid reference
  chi_mesh::Region*  aregion = this->regions.back();
  this->grid                 = aregion->GetGrid();
  size_t num_local_cells = grid->local_cells.size();

  //=================================== Set cell importance
  if (local_cell_importance_setting.empty())
  {
    local_cell_importance.clear();
    local_cell_importance.resize(num_local_cells,1.0);
  }
  else
    local_cell_importance = local_cell_importance_setting;

  //=================================== Initialize materials
  InitMaterials();

  InitParticleBatches();

  //=================================== Initialize tallies
  InitTallies();

  //=================================== Initialize Sources
  chi_log.Log(LOG_0) << "Initializing sources";
  for (auto source : sources)
  {
    source->Initialize(grid, fv, this);
    double local_weight = source->GetParallelRelativeSourceWeight();
    for (auto& val : batch_sizes_per_loc)
      val *= local_weight*chi_mpi.process_count;
  }

  //=================================== Initialize field functions
  InitFieldFunctions();

  //=================================== Init ghost ids
  InitGhostIDs();

  //=================================== Initialize data types
  BuildMPITypes();

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}