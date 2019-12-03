#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

#include "../Source/ResidualSource/mc_rmc_source.h"

extern ChiLog chi_log;
extern ChiTimer chi_program_timer;
typedef unsigned long long TULL;


extern ChiMath chi_math_handler;

//#########################################################
/**Executes the solver*/
void chi_montecarlon::Solver::ExecuteRMCUncollided()
{
  chi_log.Log(LOG_0) << "Executing Residual MonteCarlo uncollided solver";

  chi_montecarlon::Source* src = sources.back();

  if (src->type_index != SourceTypes::RESIDUAL) return;
  chi_log.Log(LOG_0) << "Mother fucker";
  uncollided_only = true;

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
      chi_montecarlon::Particle prtcl = src->CreateParticle(&rng0);

      if (prtcl.alive) nps++;

      while (prtcl.alive and !prtcl.banked) RaytraceRMC(prtcl);

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
      }
      GetOutboundBankSize();
    }


    RendesvouzTallies();
    if (make_pwld)
      RendesvouzPWLTallies();
    ComputeRelativeStdDev();



    time = chi_program_timer.GetTime()/1000.0;
    particle_rate = ((double)nps_global)*3600.0e-6/(time-start_time);
    chi_log.Log(LOG_0)
      << chi_program_timer.GetTimeString()
      << " TFC-rendesvouz: # of particles ="
      << std::setw(14)
      << nps_global
      << " avg-rate = "
      << std::setw(6) << std::setprecision(4)
      << particle_rate << " M/hr"
      << " Max Rel.Sigma = "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error
      << " "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error2
      << " "
      << std::setw(6) << std::setprecision(4) << std::scientific
      << max_relative_error3;

  }

  //Normalize tallies
  NormalizeTallies();
  if (make_pwld) NormalizePWLTallies();

  ComputePWLDTransformations();

  uncollided_only = false;

  chi_log.Log(LOG_0) << "Done executing Residual MonteCarlo uncollided solver";
}