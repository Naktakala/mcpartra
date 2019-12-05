#include "solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

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

  //============================================= If mesh global
  if (mesh_is_global)
  {
    //============================ Merge phi_local
    int num_values = phi_tally.size();
    std::vector<double> temp_phi_global(num_values,0.0);

    MPI_Allreduce(phi_tally.data(),
                  temp_phi_global.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    for (int i=0; i<num_values; i++)
      phi_global[i] += temp_phi_global[i];

    //============================ Merge the square of phi_local
    phi_global_tally_sqr.assign(num_values,0.0);

    MPI_Allreduce(phi_tally_sqr.data(),
                  phi_global_tally_sqr.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    //============================ Reset tallies
    phi_tally.assign(phi_tally.size(),0.0);
    phi_tally_sqr.assign(phi_tally_sqr.size(),0.0);
  }
    //============================================= Mesh partitioned
  else
  {
    int num_values = phi_tally.size();

    phi_global_tally_sqr.assign(num_values,0.0);

    for (int i=0; i<num_values; i++)
    {
      phi_global[i] += phi_tally[i];
      phi_global_tally_sqr[i] = phi_tally_sqr[i];
    }

    //============================ Reset tallies
    phi_tally.assign(phi_tally.size(),0.0);
    phi_tally_sqr.assign(phi_tally_sqr.size(),0.0);
  }


}

//###################################################################
/**Merges tallies from multiple locations.*/
void chi_montecarlon::Solver::RendesvouzPWLTallies()
{
  //============================================= If mesh global
  if (mesh_is_global)
  {
    //============================ Merge phi_local
    int num_values = phi_pwl_tally.size();
    std::vector<double> temp_phi_pwl_global(num_values,0.0);

    MPI_Allreduce(phi_pwl_tally.data(),
                  temp_phi_pwl_global.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    for (int i=0; i<num_values; i++)
      phi_pwl_global[i] += temp_phi_pwl_global[i];

    //============================ Merge square if phi_local
    phi_pwl_global_tally_sqr.assign(num_values,0.0);

    MPI_Allreduce(phi_pwl_tally_sqr.data(),
                  phi_pwl_global_tally_sqr.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    //============================ Reset tallies
    phi_pwl_tally.assign(phi_pwl_tally.size(),0.0);
    phi_pwl_tally_sqr.assign(phi_pwl_tally_sqr.size(),0.0);
  }
    //============================================= Mesh partitioned
  else
  {
    int num_values = phi_pwl_tally.size();

    phi_pwl_global_tally_sqr.assign(num_values,0.0);

    for (int i=0; i<num_values; i++)
    {
      phi_pwl_global[i] += phi_pwl_tally[i];
      phi_pwl_global_tally_sqr[i] = phi_pwl_tally_sqr[i];
    }

    //============================ Reset tallies
    phi_pwl_tally.assign(phi_pwl_tally.size(),0.0);
    phi_pwl_tally_sqr.assign(phi_pwl_tally_sqr.size(),0.0);
  }
}