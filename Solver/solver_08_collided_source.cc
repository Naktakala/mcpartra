#include "solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog& chi_log;

#include <chi_mpi.h>
extern ChiMPI& chi_mpi;

//###################################################################
/**Develop the collided source.*/
void chi_montecarlon::Solver::DevelopCollidedSource()
{
  size_t num_local_cells = grid->local_cell_glob_indices.size();
  cell_residual_cdf.resize(num_local_cells,0.0);

  double total_uncollided = 0.0;
  for (size_t lc=0; lc<num_local_cells; ++lc)
    total_uncollided += std::fabs(phi_uncollided_rmc[lc]);

  double running_total = 0.0;
  double volume = 0.0;
  for (size_t lc=0; lc<num_local_cells; ++lc)
  {
    running_total += std::fabs(phi_uncollided_rmc[lc]);
    cell_residual_cdf[lc] = running_total/total_uncollided;

    auto fv_view = fv_discretization->MapFeView(lc);
    volume += fv_view->volume;
  }

  MPI_Allreduce(&volume,&domain_volume,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  for (auto& val : batch_sizes_per_loc)
    val *= (volume/domain_volume)*chi_mpi.process_count;

}