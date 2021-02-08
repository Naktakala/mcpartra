#include "solver_montecarlon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/** The montecarlo method requires a particle to transfer from
 * one partition to another. This means that its local id might change. */
void chi_montecarlon::Solver::InitGhostIDs()
{
  chi_log.Log(LOG_0) << "Initializing Ghost IDs";

  //======================================== Develop list of cells requiring
  //                                         local ids
  std::vector<std::vector<uint64_t>> cells_needing_local_ids;
  cells_needing_local_ids.resize(chi_mpi.process_count);
  for (auto& cell : grid->local_cells)
    for (auto& face : cell.faces)
      if (face.has_neighbor and (not face.IsNeighborLocal(*grid)))
        cells_needing_local_ids[face.GetNeighborPartitionID(*grid)].push_back(face.neighbor_id);


  //======================================== Send counts and displacements
  std::vector<int> send_counts(chi_mpi.process_count,0);
  std::vector<int> send_displs(chi_mpi.process_count,0);
  size_t tot_send=0;
  for (int loc=0; loc<chi_mpi.process_count; ++loc)
  {
    send_counts[loc] = cells_needing_local_ids[loc].size();
    send_displs[loc] = tot_send;

    tot_send += cells_needing_local_ids[loc].size();
  }

  //======================================== Communicate receive counts
  std::vector<int> recv_counts(chi_mpi.process_count,0);
  MPI_Alltoall(send_counts.data(),1,MPI_INT,
               recv_counts.data(),1,MPI_INT,MPI_COMM_WORLD);


  //======================================== Make serial send buffer
  std::vector<uint64_t> send_buf;
  send_buf.resize(tot_send,0);
  size_t count=0;
  for (int loc=0; loc<chi_mpi.process_count; ++loc)
    for (int cell_g_id : cells_needing_local_ids[loc])
    {
      send_buf[count] = cell_g_id;
      ++count;
    }

  //======================================== Total recv size and recv displs
  int tot_recv = 0;
  std::vector<int> recv_displs(chi_mpi.process_count,0);
  for (int loc=0; loc<chi_mpi.process_count; ++loc)
  {
    recv_displs[loc] = tot_recv;
    tot_recv += recv_counts[loc];
  }

  //======================================== Communicate
  std::vector<uint64_t> recv_buf;
  recv_buf.resize(tot_recv,0);

  MPI_Alltoallv(send_buf.data(),send_counts.data(),send_displs.data(),MPI_UNSIGNED_LONG_LONG,
                recv_buf.data(),recv_counts.data(),recv_displs.data(),MPI_UNSIGNED_LONG_LONG,
                MPI_COMM_WORLD);

  //We now have all the global ids from each location
  //which needs local ids.
  std::vector<uint64_t> locl_mapping(tot_recv,0);
  auto& local_cell_glob_indices = grid->local_cell_glob_indices;
  for (int v=0; v<tot_recv; ++v)
  {
    auto location_in_native_cellstack =
      std::find(local_cell_glob_indices.begin(),
                local_cell_glob_indices.end(),
                recv_buf[v]);
    if (location_in_native_cellstack != local_cell_glob_indices.end())
      locl_mapping[v] = location_in_native_cellstack -
                        local_cell_glob_indices.begin();
    else
      chi_log.Log(LOG_ALLWARNING)
        << "chi_montecarlon::Solver::InitGhostIDs "
        << "local mapping failed.";
  }

  //======================================== Inverse notations
  auto& rsend_counts = recv_counts;
  auto& rsend_displs = recv_displs;
  auto& rrecv_counts = send_counts;
  auto& rrecv_displs = send_displs;

  std::vector<uint64_t> mappings(tot_send,-1);

  MPI_Alltoallv(locl_mapping.data(),
                rsend_counts.data(),
                rsend_displs.data(),
                MPI_UNSIGNED_LONG_LONG,
                mappings.data(),
                rrecv_counts.data(),
                rrecv_displs.data(),
                MPI_UNSIGNED_LONG_LONG,
                MPI_COMM_WORLD);

  //======================================== Create map
  for (int v=0; v<tot_send; ++v)
    cell_neighbor_nonlocal_local_id[send_buf[v]] = mappings[v];

}