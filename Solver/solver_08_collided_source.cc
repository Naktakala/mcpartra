#include "solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**Develop the collided source.*/
void chi_montecarlon::Solver::DevelopCollidedSource()
{
  size_t num_local_cells = grid->local_cell_glob_indices.size();
  cell_residual_cdf.resize(num_local_cells,0.0);

  double total_uncollided = 0.0;
  for (size_t lc=0; lc<num_local_cells; ++lc)
    total_uncollided += std::fabs(phi_uncollided_rmc[lc]);



  chi_log.Log(LOG_0) << "Collided source cdf: " << total_uncollided;
  double running_total = 0.0;
  for (size_t lc=0; lc<num_local_cells; ++lc)
  {
    running_total += std::fabs(phi_uncollided_rmc[lc]);
    cell_residual_cdf[lc] = running_total/total_uncollided;
    chi_log.Log(LOG_0) << cell_residual_cdf[lc];
  }




}