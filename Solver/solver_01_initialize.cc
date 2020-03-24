#include"solver_montecarlon.h"
#include <ChiPhysics/chi_physics.h>

#include <ChiMesh/MeshHandler/chi_meshhandler.h>

#include <ChiMath/Statistics/cdfsampler.h>

extern ChiPhysics chi_physics_handler;

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

//###################################################################
/**Initialize the solver.*/
bool chi_montecarlon::Solver::Initialize()
{
  chi_log.Log(LOG_0) << "Initializing MonteCarlo solver.";

  //=================================== Obtain grid reference
  chi_mesh::Region*  aregion = this->regions.back();
  this->grid                 = aregion->GetGrid();

  //=================================== Initialize materials
  InitMaterials();

  InitParticleBatches();

  //=================================== Initialize tallies
  InitTallies();

  //=================================== Initialize Sources
  chi_log.Log(LOG_0) << "Initializing sources";
  for (auto source : sources)
  {
    source->Initialize(grid,fv_discretization,this);
    double local_weight = source->GetParallelRelativeSourceWeight();
    for (auto& val : batch_sizes_per_loc)
      val *= local_weight*chi_mpi.process_count;
  }

  //=================================== Initialize field functions
  InitFieldFunctions();

  //======================================== Initialize data types
  BuildMPITypes();

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}