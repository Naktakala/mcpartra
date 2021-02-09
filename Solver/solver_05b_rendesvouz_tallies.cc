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

  //============================================= All grid tally rendesvouz
  for (auto& tally : grid_tally_blocks)
  {
    size_t tally_size = tally.tally_local.size();

    for (size_t i=0; i<tally_size; i++)
    {
      tally.tally_global[i]    += tally.tally_local[i];
      tally.tally_sqr_global[i] += tally.tally_sqr_local[i];
    }

    //============================ Reset tallies
    tally.tally_local.assign(tally.tally_local.size(),0.0);
    tally.tally_sqr_local.assign(tally.tally_sqr_local.size(),0.0);
  }

  //============================================= Custom tally rendesvouz
  for (auto& custom_tally : custom_tallies)
  {
    auto& tally = custom_tally.grid_tally;

    size_t tally_size = tally.tally_local.size();

    for (size_t i=0; i<tally_size; i++)
    {
      double local_val_value = tally.tally_local[i];
      double local_sqr_value = tally.tally_sqr_local[i];

      double global_val_value = 0.0;
      double global_sqr_value = 0.0;

      MPI_Allreduce(&local_val_value,  //sendbuf
                    &global_val_value, //recvbuf
                    1,                 //count
                    MPI_DOUBLE,        //datatype
                    MPI_SUM,           //operation
                    MPI_COMM_WORLD);   //communicator

      MPI_Allreduce(&local_sqr_value,  //sendbuf
                    &global_sqr_value, //recvbuf
                    1,                 //count
                    MPI_DOUBLE,        //datatype
                    MPI_SUM,           //operation
                    MPI_COMM_WORLD);   //communicator

      tally.tally_global[i]    += global_val_value;
      tally.tally_sqr_global[i] += global_sqr_value;
    }

    //============================ Reset tallies
    tally.tally_local.assign(tally.tally_local.size(),0.0);
    tally.tally_sqr_local.assign(tally.tally_sqr_local.size(),0.0);
  }
}