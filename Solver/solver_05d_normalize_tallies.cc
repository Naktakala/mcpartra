#include "solver_montecarlon.h"

//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::NormalizeTallies()
{
  if (nps_global == 0) nps_global = 1;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];
    auto cell_fv_view = fv_discretization->MapFeView(cell->cell_global_id);

    double V = cell_fv_view->volume;

    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=hi; g<=lo; g++)
    {
      int ir = lc*num_grps + g;

      phi_global[ir] *= tally_multipl_factor/nps_global/V;
      phi_global[ir] += phi_global_initial_value[ir];
//      if (not phi_uncollided_rmc.empty())
//        phi_global[ir] += phi_uncollided_rmc[ir];
    }//for g
  }//for local cell
}


//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::NormalizePWLTallies()
{
  if (nps_global == 0) nps_global = 1;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell_pwl_view   = pwl_discretization->MapFeView(cell_glob_index);
    int map             = local_cell_pwl_dof_array_address[lc];

    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=hi; g<=lo; g++)
    {
      for (int dof=0; dof<cell_pwl_view->dofs; dof++)
      {
        int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
        double V = cell_pwl_view->IntV_shapeI[dof];

        phi_pwl_global[ir] *= tally_multipl_factor/nps_global/V;

      }//for dof
    }//for g
  }//for local cell
}