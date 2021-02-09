#include "solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>
#include <ChiTimer/chi_timer.h>

#include "../Source/ResidualSource/mc_rmcA_source.h"
#include "../Source/ResidualSource/mc_rmcB_source.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiTimer chi_program_timer;

void chi_montecarlon::Solver::PrintBatchInfo(int b, double particle_rate)
{
  if (b==0 or ((b+1)%10)==0)
    chi_log.Log(LOG_0)
      << "\n"
      << "                   # of particles  prtcl-rate max-sigma  max-sigma  avg-sigma  max-fem-si \n"
      << "                                    [M/hr]    absolute   relative   absolute   absolute";
//    << "00:01:17 Batch  10       1000000      57.78   3.6438e-04 1.0356e-01 6.5527e-05";

  double correction = 1.0;
  if (typeid(*sources.back()) == typeid(chi_montecarlon::ResidualSourceA))
    correction = 1/3.0;

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Batch "
    << std::setw(3)
    << b + 1
    << std::setw(14)
    << nps_global
    << "   "
    << std::setw(8) << std::setprecision(4)
    << particle_rate
    << "   "
    << std::setw(10) << std::setprecision(4) << std::scientific
    << max_sigma*source_normalization*correction
    << " "
    << std::setw(10) << std::setprecision(4) << std::scientific
    << max_relative_sigma*source_normalization*correction
    << " "
    << std::setw(10) << std::setprecision(4) << std::scientific
    << avg_sigma*source_normalization*correction
    << " "
    << std::setw(10) << std::setprecision(4) << std::scientific
    << max_fem_sigma*source_normalization*correction;
}