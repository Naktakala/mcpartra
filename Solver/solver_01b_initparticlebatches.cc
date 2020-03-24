#include"solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

typedef unsigned long long TULL;

//###################################################################
/**Initializes batch sizes for particles to be run.*/
void chi_montecarlon::Solver::InitParticleBatches()
{
  //=================================== Process rendesvous intervals
  int num_unc_batches =
    std::ceil((double)num_uncollided_particles/tally_rendezvous_intvl);
  uncollided_batch_sizes.resize(num_unc_batches,tally_rendezvous_intvl);
  uncollided_batch_sizes[num_unc_batches-1] = num_uncollided_particles -
                                              (num_unc_batches-1)*tally_rendezvous_intvl;

  //=================================== Procces num_part per interval
  for (int b=0; b<num_unc_batches; b++)
  {
    TULL loc_num_part =
      std::ceil((double)uncollided_batch_sizes[b]/chi_mpi.process_count);

    uncollided_batch_sizes_per_loc.push_back(loc_num_part);
    if (chi_mpi.location_id == (chi_mpi.process_count-1))
      uncollided_batch_sizes_per_loc[b] = uncollided_batch_sizes[b] -
                                          chi_mpi.location_id*loc_num_part;
  }

  //=================================== Process rendesvous intervals
  int num_batches = std::ceil((double)num_particles/tally_rendezvous_intvl);
  batch_sizes.resize(num_batches,tally_rendezvous_intvl);
  batch_sizes[num_batches-1] = num_particles -
                               (num_batches-1)*tally_rendezvous_intvl;

  chi_log.Log(LOG_0) << "Number of MPI merges to be performed: " << num_batches;

  //=================================== Procces Number of particles per location
  for (int b=0; b<num_batches; b++)
  {
    TULL loc_num_part = std::ceil((double)batch_sizes[b]/chi_mpi.process_count);

    batch_sizes_per_loc.push_back(loc_num_part);
    if (chi_mpi.location_id == (chi_mpi.process_count-1))
    {
      batch_sizes_per_loc[b] = batch_sizes[b] -
                               chi_mpi.location_id*loc_num_part;
    }
  }
  chi_log.Log(LOG_0) << "Batches seperated among processes.";
  MPI_Barrier(MPI_COMM_WORLD);
}