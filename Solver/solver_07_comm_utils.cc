#include "solver_montecarlon.h"

#include <chi_mpi.h>
extern ChiMPI& chi_mpi;

#include <chi_log.h>
extern ChiLog& chi_log;

#include <sstream>

//###################################################################
/**Populates MPI datatypes*/
void mcpartra::Solver::BuildMPITypes()
{
  mcpartra::Particle::BuildMPIDatatype(mpi_prtcl_data_type);
}

//###################################################################
/**Determine Outbound bank size*/
void mcpartra::Solver::GetOutboundBankSize()
{
  unsigned int local_size = outbound_particle_bank.size();

  MPI_Allreduce(&local_size,&total_outbound_bank_size,1,
                MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
}

//###################################################################
/**Builds the outgoing/incoming structures.*/
void mcpartra::Solver::
  ReceiveIncomingParticles(std::vector<Particle>& inbound_particles)
{
  //=================================== Segregate particles acc. to location
  std::vector<std::vector<Particle>> locI_outbound(chi_mpi.process_count);

  for (auto& outb_prtcl : outbound_particle_bank)
  {
    int loc_i = grid->cells[outb_prtcl.cur_cell_global_id].partition_id;

    locI_outbound[loc_i].push_back(outb_prtcl);
  }
  outbound_particle_bank.clear();

  //=================================== Determine send counts
  std::vector<int> send_counts(chi_mpi.process_count,0);
  int total_outbound = 0;
  for (int i=0; i<chi_mpi.process_count; i++)
  {
    send_counts[i] = locI_outbound[i].size();
    total_outbound += send_counts[i];
  }

  //=================================== Determine displacements and
  //                                    linearize storage
  std::vector<int> send_displacements(chi_mpi.process_count,0);
  std::vector<Particle> linear_outbound;
  linear_outbound.reserve(total_outbound);
  for (int i=0; i<chi_mpi.process_count; i++)
  {
    send_displacements[i] = linear_outbound.size();
    for (auto& outb_prtcl : locI_outbound[i])
      linear_outbound.push_back(outb_prtcl);//std::move(outb_prtcl));
  }

//  if (not linear_outbound.empty())
//  {
//    auto& prtl0 = linear_outbound[0];
//    chi_log.Log(LOG_ALL)
//      << "Send "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }
//  if (not linear_outbound.empty())
//  {
//    auto& prtl0 = linear_outbound[1];
//    chi_log.Log(LOG_ALL)
//      << "Send "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }
//  if (not linear_outbound.empty())
//  {
//    auto& prtl0 = linear_outbound[2];
//    chi_log.Log(LOG_ALL)
//      << "Send "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }
//  chi_log.Log(LOG_ALL) << "Communicating";

  //=================================== Collect recv counts
  std::vector<int> recv_counts(chi_mpi.process_count,0);
  MPI_Alltoall(send_counts.data(),1,MPI_INT,
               recv_counts.data(),1,MPI_INT,
               MPI_COMM_WORLD);

  //=================================== Determine recv displacements
  std::vector<int> recv_displacements(chi_mpi.process_count,0);
  int total_inbound = 0;
  for (int i=0; i<chi_mpi.process_count; i++)
  {
    recv_displacements[i] = total_inbound;
    total_inbound += recv_counts[i];

//    chi_log.Log(LOG_ALL)
//      << "Recv from " << i
//      << " " << recv_counts[i]
//      << " " << recv_displacements[i];
  }

//  std::vector<Particle> linear_inbound(total_inbound,Particle());
  inbound_particles.resize(total_inbound,Particle());

  int error_code = MPI_Alltoallv(
    linear_outbound.data(),
    send_counts.data(),
    send_displacements.data(),mpi_prtcl_data_type,
    inbound_particles.data(),
    recv_counts.data(),
    recv_displacements.data(),mpi_prtcl_data_type,
    MPI_COMM_WORLD);

  //======================================== Add source particles
//  for (auto& prtcl : particle_source_bank)
//    inbound_particles.push_back(prtcl);
//
//  particle_source_bank.clear();

//  if (error_code != MPI_SUCCESS)
//  {
//    std::stringstream err_stream;
//    err_stream << "All to all error.";
//    char error_string[BUFSIZ];
//    int length_of_error_string, error_class;
//    MPI_Error_class(error_code, &error_class);
//    MPI_Error_string(error_class, error_string, &length_of_error_string);
//    err_stream << error_string << "\n";
//    MPI_Error_string(error_code, error_string, &length_of_error_string);
//    err_stream << error_string << "\n";
//    chi_log.Log(LOG_ALLWARNING) << err_stream.str();
//  }

//  if (total_inbound>0)
//  {
//    auto& prtl0 = inbound_particles[0];
//    chi_log.Log(LOG_ALL)
//      << "Recv "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }
//  if (total_inbound>0)
//  {
//    auto& prtl0 = inbound_particles[1];
//    chi_log.Log(LOG_ALL)
//      << "Recv "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }
//  if (total_inbound>0)
//  {
//    auto& prtl0 = inbound_particles[2];
//    chi_log.Log(LOG_ALL)
//      << "Recv "
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.egrp << " "
//      << prtl0.cur_cell_global_id << " "
//      << prtl0.pre_cell_global_id << " "
//      << prtl0.cur_cell_local_id << " "
//      << prtl0.pre_cell_local_id << " "
//      << prtl0.alive << " "
//      << prtl0.banked;
//  }

}