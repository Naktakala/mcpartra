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

  //TODO: Begin Remove test code
  size_t num_stored_particles = options.num_particles;
  std::vector<mcpartra::Particle> source_particles;
  source_particles.reserve(num_stored_particles);

  chi_math::RandomNumberGenerator rng1;

  chi_mesh::Vector3 avg_pos, max_pos, min_pos(100,100,100);
  for (unsigned long long b : batch_sizes_per_loc)
  {
    for (uint64_t pi=0; pi<b; ++pi)
    {
      mcpartra::Particle prtcl = SampleSources(rng1);
      avg_pos += prtcl.dir;

      max_pos.x = std::max(max_pos.x,prtcl.pos.x);
      max_pos.y = std::max(max_pos.y,prtcl.pos.y);
      max_pos.z = std::max(max_pos.z,prtcl.pos.z);

      min_pos.x = std::min(min_pos.x,prtcl.pos.x);
      min_pos.y = std::min(min_pos.y,prtcl.pos.y);
      min_pos.z = std::min(min_pos.z,prtcl.pos.z);

      source_particles.push_back(prtcl);
//      chi_log.Log() << prtcl.dir.PrintS();
    }//for pi in batch
  }//for batch
  size_t k=0;
  avg_pos /= source_particles.size();
  chi_log.Log(LOG_0) << "Average source position: " << avg_pos.PrintS();
  chi_log.Log(LOG_0) << "Max position: " << max_pos.PrintS();
  chi_log.Log(LOG_0) << "Min position: " << min_pos.PrintS();
  //TODO: End Remove test code

//  if (TextName() == "FMCParTra")
//    WriteParticlesToFile("ZParticles.data", source_particles);
//  if (TextName() == "RMCParTra")
//    source_particles = ReadParticlesFromFile("ZParticles.data");

  //============================================= Start batch processing
  double start_time = chi_program_timer.GetTime()/1000.0;
  size_t start_nps_global = nps_global;

  for (size_t b=0; b<batch_sizes_per_loc.size(); b++)
  {
    nps = 0;
    for (uint64_t pi=0; pi<batch_sizes_per_loc[b]; ++pi)
    {
      mcpartra::Particle prtcl = SampleSources(rng0);
//      mcpartra::Particle prtcl = source_particles[k++];

      if (prtcl.alive) nps++;

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

    PrintBatchInfo(b,particle_rate);
  }//for batch

  //============================================= Post processing
  if (options.write_run_tape) WriteRunTape(options.run_tape_base_name);

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
    {
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
    }//for m
    outstr << "\n";

    chi_log.Log() << outstr.str();
  }

  //============================================= Accumulated source importances
  double accumulated_src_importances_global;
  MPI_Allreduce(&accumulated_src_importances,        //sendbuf
                &accumulated_src_importances_global, //recvbuf
                1,                                   //count
                MPI_DOUBLE,                          //datatype
                MPI_SUM,                             //operation
                MPI_COMM_WORLD);                     //communicator
  chi_log.Log() << "MCParTra: average source particle importance = "
                << std::scientific
                << accumulated_src_importances_global/
                   static_cast<double>(nps_global);

  chi_log.Log(LOG_0) << "Done executing Montecarlo solver";
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