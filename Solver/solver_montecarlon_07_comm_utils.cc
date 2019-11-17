#include "solver_montecarlon.h"

#include <chi_mpi.h>
extern ChiMPI chi_mpi;

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**Populates MPI datatypes*/
void chi_montecarlon::Solver::BuildMPITypes()
{
  int block_lengths[] = {3,3,1,2,2};
  MPI_Aint block_disp[] = {
    offsetof(Particle,pos),
    offsetof(Particle,dir),
    offsetof(Particle,w),
    offsetof(Particle,egrp),
    offsetof(Particle,alive)};

  MPI_Datatype block_types[] = {
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_DOUBLE,
    MPI_INT,
    MPIU_BOOL
  };

  MPI_Type_create_struct(
    5,              // Number of blocks
    block_lengths,  // Block lengths
    block_disp,     // Block displacements
    block_types,    // Block types
    &mpi_prtcl_data_type);

  MPI_Type_commit(&mpi_prtcl_data_type);
}

//###################################################################
/**Determine Outbound bank size*/
void chi_montecarlon::Solver::GetOutboundBankSize()
{
  unsigned int local_size = outbound_particle_bank.size();

  MPI_Allreduce(&local_size,&total_outbound_bank_size,1,
                MPI_UNSIGNED,MPI_SUM,MPI_COMM_WORLD);
}

//###################################################################
/**Builds the outgoing/incoming structures.*/
void chi_montecarlon::Solver::
  ReceiveIncomingParticles(std::vector<Particle>& inbound_particles)
{
  //=================================== Segregate particles acc. to location
  std::vector<std::vector<Particle>> locI_outbound(chi_mpi.process_count);

  for (auto& outb_prtcl : outbound_particle_bank)
  {
    int loc_i = grid->cells[outb_prtcl.cur_cell_ind]->partition_id;

    locI_outbound[loc_i].push_back(outb_prtcl); //std::move(outb_prtcl));
  }
  outbound_particle_bank.clear();

  //=================================== Determine send counts
  std::vector<int> send_counts(chi_mpi.process_count,0);
  int total_outbound = 0;
  for (int i=0; i<chi_mpi.process_count; i++)
  {
    send_counts[i] = locI_outbound[i].size();
    total_outbound += send_counts[i];
//    chi_log.Log(LOG_ALL)
//      << "Send to " << i
//      << " " << send_counts[i];
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
//      << prtl0.pos.PrintS() << " "
//      << prtl0.dir.PrintS() << " "
//      << prtl0.w << " "
//      << prtl0.cur_cell_ind << " "
//      << prtl0.egrp;
//  }

}