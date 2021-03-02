#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

extern ChiLog& chi_log;
extern ChiTimer chi_program_timer;
typedef unsigned long long TULL;

#include "../Source/ResidualSource/mc_rmcA_source.h"

//#########################################################
/**Executes the solver*/
void chi_montecarlon::Solver::Execute()
{
  chi_log.Log(LOG_0) << "Executing Montecarlo solver";

  chi_montecarlon::Source* src = sources.back();

  std::vector<Particle> inbound_particles;

  nps_global = 0;
  double start_time = chi_program_timer.GetTime()/1000.0;
  double time = 0.0;
  double particle_rate = 0.0;
  for (size_t b=0; b<batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    current_batch = b;

    for (TULL pi=0; pi<batch_sizes_per_loc[b]; pi++)
    {
//      chi_log.Log() << "src";
      chi_montecarlon::Particle prtcl = src->CreateParticle(&rng0);

      if (prtcl.alive) nps++;
//      prtcl.alive = false;

      while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

      while (not particle_source_bank.empty())
      {
        prtcl = particle_source_bank.back();
        particle_source_bank.pop_back();

        prtcl.alive = true;
        prtcl.banked = false;

        while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);
      }
    }//for pi in batch

    MPI_Barrier(MPI_COMM_WORLD);

    GetOutboundBankSize();
    while (total_outbound_bank_size>0)
    {
      ReceiveIncomingParticles(inbound_particles);
      for (auto& prtcl : inbound_particles)
      {
        prtcl.alive = true;
        prtcl.banked = false;
        while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);

        while (not particle_source_bank.empty())
        {
          prtcl = particle_source_bank.back();
          particle_source_bank.pop_back();

          prtcl.alive = true;
          prtcl.banked = false;

          while (prtcl.alive and !prtcl.banked) Raytrace(prtcl);
        }
      }
      GetOutboundBankSize();
    }


    RendesvouzTallies();
    ComputeUncertainty();

    time = chi_program_timer.GetTime()/1000.0;
    particle_rate = ((double)nps_global)*3600.0e-6/(time-start_time);

    PrintBatchInfo(b,particle_rate);
  }//for batch

  //Normalize tallies
  NormalizeTallies();

  ComputePWLDTransformations();

  double correction=1.0;
  {
    if (typeid(*sources.back()) == typeid(chi_montecarlon::ResidualSourceA))
      correction = 1/3.0;
  }

  //============================================= Print custom tallies TFC
  int cust_counter = -1;
  for (auto& custom_tally : custom_tallies)
  {
    std::stringstream outstr;
    outstr << "Custom tally " << ++cust_counter << " TFC:\n";

    int comp_counter = -1;
    for (int m=0; m<num_moms; ++m)
      for (int g=0; g<num_grps; ++g)
      {
        outstr << "Component " << ++comp_counter << ":\n";

        auto ir = uk_man_fv.MapUnknown(m, g);

        for (auto& tfc : custom_tally.tally_fluctation_chart)
        {
          outstr
            << std::setw(10) << std::setprecision(4) << std::scientific
            << tfc.average[ir]*correction
            << " "
            << std::setw(10) << std::setprecision(4) << std::scientific
            << tfc.sigma[ir]*correction << "\n";
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