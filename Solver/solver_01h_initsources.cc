#include "solver_montecarlon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Initialize Sources.*/
void mcpartra::Solver::InitSources()
{
  chi_log.Log() << "MCParTra: Initializing sources";

  for (auto& source : sources)
  {
    source->Initialize(grid, fv);

    total_globl_source_rate += source->GlobalSourceRate();
    total_local_source_rate += source->LocalSourceRate();
  }
  source_normalization = total_globl_source_rate;

  char buffer[100];
  snprintf(buffer,100,"%g",total_globl_source_rate);
  std::string stotal_globl_src_rate(buffer);

  chi_log.Log() << "  Total amount of sources     = " << sources.size();
  chi_log.Log() << "  Total source rate           = " << stotal_globl_src_rate;

  //============================================= Build local cdf
  size_t num_sources = sources.size();
  local_source_cdf.assign(num_sources, 0.0);
  double running_total = 0.0;
  for (int s=0; s<num_sources; ++s)
  {
    running_total += sources[s]->LocalSourceRate();
    local_source_cdf[s] = running_total/total_local_source_rate;
  }

}