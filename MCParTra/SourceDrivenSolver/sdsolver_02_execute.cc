#include "sdsolver.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"

extern ChiLog& chi_log;
extern ChiTimer chi_program_timer;

#include "Sources/ResidualSource/mc_rmcA_source.h"

//###################################################################
/**Executes the solver*/
void mcpartra::SourceDrivenSolver::Execute()
{
  chi_log.Log(LOG_0) << "\nExecuting MCParTra solver\n";

  //============================================= Start batch processing
  double start_time = chi_program_timer.GetTime()/1000.0;
  size_t start_nps_global = nps_global;

  double avg_particle_rate = 0.0;
  for (size_t b=0; b<batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    for (uint64_t pi=0; pi<batch_sizes_per_loc[b]; ++pi)
    {
      mcpartra::Particle prtcl = SampleSources(rng0);

      if (prtcl.alive) nps++;
      if (options.no_transport) prtcl.Kill();

      while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

      DepleteSourceParticleBank();
    }//for pi in batch

    //================================= Communicate inbound/outbound particles
    //                                  and deplete banks
    CommAndSimOutboundInboundParticles();

    //================================= Rendesvouz and housekeeping
    RendesvouzTallies();
    ComputeUncertainty();

    double time = chi_program_timer.GetTime()/1000.0;
    double particle_rate = ((double)(nps_global-start_nps_global)*3.6e-3/
                           (time-start_time));

    avg_particle_rate += particle_rate;
    PrintBatchInfo(b,particle_rate);
  }//for batch
  avg_particle_rate /= batch_sizes_per_loc.size();

//  exit(EXIT_SUCCESS);

  //============================================= Post processing
  if (options.write_run_tape) WriteRunTape(options.run_tape_base_name);

  chi_log.Log() << "\nMCParTra: Ray tracing complete.\n\n";

  chi_log.Log() << "MCParTra: Number of lost particles = "
                << lost_particles.size();
  if (chi_log.GetVerbosity() > LOG_0VERBOSE_2)
  {
    std::stringstream lost_particle_log;
    for (const auto& lost_particle : lost_particles)
      lost_particle_log << lost_particle << "\n";
    chi_log.Log(LOG_ALLVERBOSE_2)
    << "MCParTra: Lost particle log:" << lost_particle_log.str();
  }

  NormalizeTallies();
  ComputePWLDTransformations();

  chi_log.Log() << "MCParTra: Average particle simulation rate: "
                << avg_particle_rate << " [M/hr].";
  PrintCustomTallies();

  chi_log.Log(LOG_0) << "\nDone executing Montecarlo solver\n\n";
}

//###################################################################
/**Short method to deplete bank.*/
void mcpartra::SourceDrivenSolver::DepleteSourceParticleBank()
{
  mcpartra::Particle prtcl;
  while (not particle_source_bank.empty())
  {
    prtcl = particle_source_bank.back();
    particle_source_bank.pop_back();

    prtcl.alive = true;
    prtcl.banked = false;

    while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);
  }
}

//###################################################################
/**Communicate outbound and inbound particles, then
 * deplete the inbound banks.*/
void mcpartra::SourceDrivenSolver::CommAndSimOutboundInboundParticles()
{
  std::vector<Particle> inbound_particle_bank;

  GetOutboundBankSize();
  while (total_outbound_bank_size>0)
  {
    ReceiveIncomingParticles(inbound_particle_bank);
    for (auto& prtcl : inbound_particle_bank)
    {
      prtcl.alive = true;
      prtcl.banked = false;
      while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

      DepleteSourceParticleBank();
    }
    GetOutboundBankSize();
  }
}