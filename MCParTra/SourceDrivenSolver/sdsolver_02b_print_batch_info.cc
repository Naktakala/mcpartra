#include "sdsolver.h"

#include "chi_log.h"
#include "chi_mpi.h"
#include "ChiTimer/chi_timer.h"

#include "Sources/ResidualSource/mc_rmcA_source.h"

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;
extern ChiTimer chi_program_timer;

void mcpartra::SourceDrivenSolver::PrintBatchInfo(size_t b, double particle_rate)
{
  if (b==0 or (b%10)==0)
    chi_log.Log(LOG_0)
      << "\n"
      << "                   # of particles  prtcl-rate\n"
      << "                                    [M/hr]    ";
//    << "00:01:17 Batch  10       1000000      57.78

  chi_log.Log(LOG_0)
    << chi_program_timer.GetTimeString()
    << " Batch "
    << std::setw(3)
    << b + 1
    << std::setw(14)
    << nps_global
    << "   "
    << std::setw(8) << std::setprecision(4)
    << particle_rate;
}

void mcpartra::SourceDrivenSolver::PrintCustomTallies()
{
  double correction=1.0;
  {
    if (sources.back()->Type() == RESIDUAL_TYPE_A)
      correction = 1/3.0;
  }

  int cust_counter = 0;
  for (auto& tally : custom_tallies)
  {
    std::stringstream outstr;

    //=========================================== Print TFC
    if (options.print_TFC)
    {
      outstr << "Custom tally " << cust_counter++ << " TFC:\n";
      for (size_t m=0; m < num_moments; ++m)
      {
        for (size_t g=0; g < num_groups; ++g)
        {
          outstr << "m=" << m << " g=" << " :\n";
          auto dof_map = uk_man_fv.MapUnknown(m, g);

          for (auto& tfc : tally.tally_fluctation_chart)
          {
            tfc.sigma[dof_map] *= correction;
            const double avg = tfc.average[dof_map];
            const double stddev = tfc.sigma[dof_map];
            const double relative_stddev = (std::isinf(stddev/avg) or
                                            std::isnan(stddev/avg))? 0.0 :
                                           stddev/avg;

            outstr
              << std::setw(10) << std::setprecision(4) << std::scientific
              << avg << " "
              << std::setw(10) << std::setprecision(4) << std::scientific
              << stddev << " "
              << std::setw(10) << std::setprecision(4) << std::scientific
              << relative_stddev << "\n";
          }
        }//for g
      }//for m
      outstr << "\n";
    }//if options.print_TFC

    //=========================================== Print final result
    outstr << "Custom tally " << cust_counter++ << " Final result:\n";
    for (size_t m=0; m < num_moments; ++m)
    {
      for (size_t g=0; g < num_groups; ++g)
      {
        outstr << "m=" << m << " g=" << g << " : ";
        auto dof_map = uk_man_fv.MapUnknown(m, g);


        auto& tfc = tally.tally_fluctation_chart.back();
        {
          const double avg = tfc.average[dof_map];
          const double stddev = tfc.sigma[dof_map];
          const double relative_stddev = (std::isinf(stddev/avg) or
                                          std::isnan(stddev/avg))? 0.0 :
                                         stddev/avg;

          if (not options.print_TFC) tfc.sigma[dof_map] *= correction;
          outstr
            << std::setw(10) << std::setprecision(4) << std::scientific
            << avg << " "
            << std::setw(10) << std::setprecision(4) << std::scientific
            << stddev << " "
            << std::setw(10) << std::setprecision(4) << std::scientific
            << relative_stddev << "\n";
        }
      }//for g
    }//for m
    outstr << "\n";

    chi_log.Log() << outstr.str();
  }
}