#include "solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

typedef unsigned long long TULL;

//###################################################################
/**Merges tallies from multiple locations.*/
void chi_montecarlon::Solver::RendesvouzTallies()
{
  TULL temp_nps_global = 0;
  MPI_Allreduce(&nps,&temp_nps_global,
                1,MPI_UNSIGNED_LONG_LONG,
                MPI_SUM,MPI_COMM_WORLD);

  nps_global += temp_nps_global;

  //============================================= All tally rendesvouz
  for (auto& tally : grid_tally_blocks)
  {
    int tally_size = tally.tally_local.size();

    tally.tally_sqr_global.assign(tally_size,0.0);

    for (int i=0; i<tally_size; i++)
    {
      tally.tally_global[i]    += tally.tally_local[i];
      tally.tally_sqr_global[i] = tally.tally_sqr_local[i];
    }

    //============================ Reset tallies
    tally.tally_local.assign(tally.tally_local.size(),0.0);
    tally.tally_sqr_local.assign(tally.tally_sqr_local.size(),0.0);
  }
}