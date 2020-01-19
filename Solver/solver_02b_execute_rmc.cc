#include "solver_montecarlon.h"

#include <chi_log.h>
#include <ChiTimer/chi_timer.h>

#include "../Source/ResidualSource/mc_rmc2_source.h"

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

  auto rsrc = (chi_montecarlon::ResidualSource2*)src;
  uncollided_only = true;

  std::vector<Particle> inbound_particles;

  nps_global = 0;
  double start_time = chi_program_timer.GetTime()/1000.0;
  double time = 0.0;
  double particle_rate = 0.0;
  for (size_t b=0; b<uncollided_batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    current_batch = b;
    for (TULL pi=0; pi<uncollided_batch_sizes_per_loc[b]; pi++)
    {
      chi_montecarlon::Particle prtcl = rsrc->CreateBndryParticle(&rng0);

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
        while (prtcl.alive and !prtcl.banked) RaytraceRMC(prtcl);
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

  size_t avg_tally_size = phi_global.size();
  size_t pwl_tally_size = phi_pwl_tally_contrib.size();

  //======================================== Copy avg tally to uncollided
  phi_uncollided_rmc.clear();
  phi_uncollided_rmc.reserve(avg_tally_size);
  std::copy(phi_global.begin(),
            phi_global.end(),
            std::back_inserter(phi_uncollided_rmc));

  //======================================== Copy pwl tally to uncollided
  phi_pwl_uncollided_rmc.clear();
  phi_pwl_uncollided_rmc.reserve(pwl_tally_size);
  std::copy(phi_pwl_global.begin(),
            phi_pwl_global.end(),
            std::back_inserter(phi_pwl_uncollided_rmc));

  //======================================== Clear avg tallies
  phi_tally_contrib.clear();
  phi_tally.clear();
  phi_tally_sqr.clear();

  phi_global_initial_value.clear();
  phi_global.clear();
  phi_global_tally_sqr.clear();

  phi_local_relsigma.clear();

  //======================================== Clear pwl tallies
  phi_pwl_tally_contrib.clear();
  phi_pwl_tally.clear();
  phi_pwl_tally_sqr.clear();

  phi_pwl_global.clear();
  phi_pwl_global_tally_sqr.clear();

  phi_pwl_local_relsigma.clear();

  //======================================== Reset avg tallies
  phi_tally_contrib.resize(avg_tally_size,0.0);
  phi_tally.resize(avg_tally_size,0.0);
  phi_tally_sqr.resize(avg_tally_size,0.0);

  phi_global_initial_value.resize(avg_tally_size,0.0);
  phi_global.resize(avg_tally_size,0.0);
  phi_global_tally_sqr.resize(avg_tally_size,0.0);

  phi_local_relsigma.resize(avg_tally_size,0.0);

  //======================================== Reset PWL tallies
  phi_pwl_tally_contrib.resize(pwl_tally_size,0.0);
  phi_pwl_tally.resize(pwl_tally_size,0.0);
  phi_pwl_tally_sqr.resize(pwl_tally_size,0.0);

  phi_pwl_global.resize(pwl_tally_size,0.0);
  phi_pwl_global_tally_sqr.resize(pwl_tally_size,0.0);

  phi_pwl_local_relsigma.resize(pwl_tally_size,0.0);

  //======================================== Developing the collided source
  DevelopCollidedSource();


  uncollided_only = false;

  chi_log.Log(LOG_0) << "Done executing Residual MonteCarlo uncollided solver";
}