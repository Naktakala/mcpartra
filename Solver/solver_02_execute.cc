#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiLog& chi_log;
extern ChiTimer chi_program_timer;

#include "../Source/ResidualSource/mc_rmcA_source.h"

//#########################################################
/**Executes the solver*/
void mcpartra::Solver::Execute()
{
  chi_log.Log(LOG_0) << "\nExecuting MCParTra solver\n";

  mcpartra::SourceBase* src = sources.back();

  //#######################################################
  /**Short lambda to deplete bank.*/
  auto DepleteBank = [this]()
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
  };

  //============================================= Start batch processing
  std::vector<Particle> inbound_particles;
  double start_time = chi_program_timer.GetTime()/1000.0;
  size_t start_nps_global = nps_global;
  for (size_t b=0; b<batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    for (uint64_t pi=0; pi<batch_sizes_per_loc[b]; ++pi)
    {
      mcpartra::Particle prtcl = SampleSources(rng0);

      if (prtcl.alive) nps++;

      while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

      DepleteBank();
    }//for pi in batch

    MPI_Barrier(MPI_COMM_WORLD);

    //================================= Communicate and deplete banks
    GetOutboundBankSize();
    while (total_outbound_bank_size>0)
    {
      ReceiveIncomingParticles(inbound_particles);
      for (auto& prtcl : inbound_particles)
      {
        prtcl.alive = true;
        prtcl.banked = false;
        while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

        DepleteBank();
      }
      GetOutboundBankSize();
    }

    //================================= Rendesvouz and housekeeping
    RendesvouzTallies();
    ComputeUncertainty();

    double time = chi_program_timer.GetTime()/1000.0;
    double particle_rate = ((double)(nps_global-start_nps_global)*3.6e-3/
                           (time-start_time));

    PrintBatchInfo(b,particle_rate);
  }//for batch

  //============================================= Post processing
  if (options.write_run_tape) WriteRunTape(options.run_tape_base_name);
  NormalizeTallies();
  ComputePWLDTransformations();

  double correction=1.0;
  {
    if (sources.back()->Type() == RESIDUAL_TYPE_A)
      correction = 1/3.0;
  }

  //============================================= Print custom tallies TFC
  int cust_counter = -1;
  for (auto& custom_tally : custom_tallies)
  {
    std::stringstream outstr;
    outstr << "Custom tally " << ++cust_counter << " TFC:\n";

    int comp_counter = -1;
    for (size_t m=0; m < num_moments; ++m)
      for (size_t g=0; g < num_groups; ++g)
      {
        outstr << "Component " << ++comp_counter << ":\n";

        auto dof_map = uk_man_fv.MapUnknown(m, g);

        for (auto& tfc : custom_tally.tally_fluctation_chart)
        {
          outstr
            << std::setw(10) << std::setprecision(4) << std::scientific
            << tfc.average[dof_map] * correction
            << " "
            << std::setw(10) << std::setprecision(4) << std::scientific
            << tfc.sigma[dof_map] * correction << "\n";
        }
      }//for g
    outstr << "\n";

    chi_log.Log() << outstr.str();
  }

  if (src->CheckForReExecution())
    Execute();
  else
    chi_log.Log(LOG_0) << "Done executing Montecarlo solver";
}