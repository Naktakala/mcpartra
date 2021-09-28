#include"sdsolver.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include "ChiConsole/chi_console.h"
extern ChiConsole&  chi_console;

//###################################################################
/**Initialize the solver.*/
void mcpartra::SourceDrivenSolver::Initialize()
{
  chi_log.Log(LOG_0) << "\nInitializing MCParTra solver.\n\n";

  //=================================== Obtain grid reference
  grid = chi_mesh::GetCurrentHandler()->GetGrid();

  InitRaytracing();
  InitMaterials();
  InitMomentIndices();
  InitTallies();

  {
    double process_memory = chi_console.GetMemoryUsageInMB();
    double program_memory = 0.0;
    MPI_Allreduce(&process_memory, &program_memory, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double peak_process_memory = 0.0;
    MPI_Allreduce(&process_memory, &peak_process_memory, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    chi_log.Log(LOG_0) << "\nMCParTra: Solver initialized1. \n"
                       << "  Peak process memory  "
                       << std::fixed << std::setprecision(1)
                       << process_memory << "Mb\n"
                       << "  Program total memory "
                       << std::fixed << std::setprecision(1)
                       << program_memory << "Mb\n";
  }

  InitFieldFunctions();
  InitCellImportances();

  {
    double process_memory = chi_console.GetMemoryUsageInMB();
    double program_memory = 0.0;
    MPI_Allreduce(&process_memory, &program_memory, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double peak_process_memory = 0.0;
    MPI_Allreduce(&process_memory, &peak_process_memory, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    chi_log.Log(LOG_0) << "\nMCParTra: Solver initialized2. \n"
                       << "  Peak process memory  "
                       << std::fixed << std::setprecision(1)
                       << process_memory << "Mb\n"
                       << "  Program total memory "
                       << std::fixed << std::setprecision(1)
                       << program_memory << "Mb\n";
  }

  InitCellGeometryData();
  InitSources();

  {
    double process_memory = chi_console.GetMemoryUsageInMB();
    double program_memory = 0.0;
    MPI_Allreduce(&process_memory, &program_memory, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double peak_process_memory = 0.0;
    MPI_Allreduce(&process_memory, &peak_process_memory, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    chi_log.Log(LOG_0) << "\nMCParTra: Solver initialized3. \n"
                       << "  Peak process memory  "
                       << std::fixed << std::setprecision(1)
                       << process_memory << "Mb\n"
                       << "  Program total memory "
                       << std::fixed << std::setprecision(1)
                       << program_memory << "Mb\n";
  }

  InitParticleBatches();


  BuildMPITypes();

  {
    double process_memory = chi_console.GetMemoryUsageInMB();
    double program_memory = 0.0;
    MPI_Allreduce(&process_memory, &program_memory, 1, MPI_DOUBLE,
                  MPI_SUM, MPI_COMM_WORLD);

    double peak_process_memory = 0.0;
    MPI_Allreduce(&process_memory, &peak_process_memory, 1, MPI_DOUBLE,
                  MPI_MAX, MPI_COMM_WORLD);

    chi_log.Log(LOG_0) << "\nMCParTra: Solver initialized. \n"
                       << "  Peak process memory  "
                       << std::fixed << std::setprecision(1)
                       << process_memory << "Mb\n"
                       << "  Program total memory "
                       << std::fixed << std::setprecision(1)
                       << program_memory << "Mb\n";
  }

  MPI_Barrier(MPI_COMM_WORLD);
}