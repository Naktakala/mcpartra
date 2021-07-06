#include"solver_montecarlon.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

//###################################################################
/**Initialize the solver.*/
bool mcpartra::Solver::Initialize()
{
  chi_log.Log(LOG_0) << "\nInitializing MCParTra solver.\n\n";

  //=================================== Obtain grid reference
  grid = chi_mesh::GetCurrentHandler()->GetGrid();

  default_raytracer = std::make_shared<chi_mesh::RayTracer>(*grid,
                                                            1.0e-8,
                                                            1.0e-10,
                                                            1.0e5,
                                                            false);

  InitRaytracing();
  InitMaterials();
  InitMomentIndices();
  InitTallies();
  InitFieldFunctions();
  InitGhostIDs();
  InitCellImportances();
  InitSources();
  InitParticleBatches();


  BuildMPITypes();

  MPI_Barrier(MPI_COMM_WORLD);
  return true;
}