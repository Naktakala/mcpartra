#include"solver_montecarlon.h"

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog& chi_log;
extern ChiMPI& chi_mpi;

typedef unsigned long long TULL;

//###################################################################
/**Initializes batch sizes for particles to be run.*/
void mcpartra::Solver::InitParticleBatches()
{
  chi_log.Log() << "MCParTra: Initializing particle batches.";

  //################################### Lambda to create comma seperated string
  auto MakeCommaString = [](uint64_t integer)
  {
    auto s = std::to_string(integer);
    int n = static_cast<int>(s.length()) - 3;
    while (n > 0) {
      s.insert(n, ",");
      n -= 3;
    }
    return s;
  };

  //################################### Lambda to create scientific string
  auto MakeScientific10 = [](uint64_t integer)
  {
    std::stringstream ss;
    ss << std::scientific << std::setprecision(3) << static_cast<double>(integer);

    return ss.str();
  };

  //=================================== Process rendesvous intervals
  int num_batches = std::ceil(static_cast<double>(options.num_particles)/
                              static_cast<double>(options.tally_rendezvous_intvl));
  batch_sizes.resize(num_batches,options.tally_rendezvous_intvl);
  batch_sizes[num_batches-1] = options.num_particles -
                               (num_batches-1)*options.tally_rendezvous_intvl;

  chi_log.Log() << "  Number of Prtcls to run        = "
                << MakeScientific10(options.num_particles)
                << " or "
                << MakeCommaString(options.num_particles);
  chi_log.Log() << "  Tally Rendesvouz Interval      = "
                << MakeScientific10(options.tally_rendezvous_intvl)
                << " or "
                << MakeCommaString(options.tally_rendezvous_intvl);
  chi_log.Log() << "  Num of batches to rendesvouz   = " << num_batches;
  double local_source_fraction = total_local_source_rate/
                                 static_cast<double>(total_globl_source_rate);
  chi_log.Log() << "  Local source fraction          = " << local_source_fraction;

  //=================================== Procces Number of particles per location
  batch_sizes_per_loc.clear();
  batch_sizes_per_loc.resize(num_batches,0.0);
  for (int b=0; b<num_batches; b++)
  {
    batch_sizes_per_loc[b] = static_cast<uint64_t>(
      static_cast<double>(batch_sizes[b])*
      total_local_source_rate/total_globl_source_rate);
  }
  chi_log.Log(LOG_0) << "Batches seperated among processes.";
  MPI_Barrier(MPI_COMM_WORLD);

  //=================================== Process rendesvous intervals
  int num_unc_batches =
    std::ceil(static_cast<double>(options.num_uncollided_particles)/
              static_cast<double>(options.tally_rendezvous_intvl));
  uncollided_batch_sizes.resize(num_unc_batches,options.tally_rendezvous_intvl);
  uncollided_batch_sizes[num_unc_batches-1] =
    options.num_uncollided_particles -
    (num_unc_batches-1)*options.tally_rendezvous_intvl;

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
}