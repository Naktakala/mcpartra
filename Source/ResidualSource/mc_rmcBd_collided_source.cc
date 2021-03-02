#include "mc_rmcB_source.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics& chi_physics_handler;

//###################################################################
/**Develops a collided source.*/
void chi_montecarlon::ResidualSourceB::
  DevelopRMCCollidedSource()
{
  size_t num_local_cells = grid->local_cell_glob_indices.size();
  const int G = ref_solver->num_grps;

  //============================================= Initialize sizes and integrals
  ref_solver->cdf_phi_unc_group.assign(G, 0.0);
  ref_solver->cdf_phi_unc_group_cell.resize(G);
  for (int g=0; g<G; ++g)
    ref_solver->cdf_phi_unc_group_cell[g].assign(num_local_cells, 0.0);
  ref_solver->IntVSumG_phi_unc = 0.0;
  ref_solver->IntV_phi_unc_g.assign(G,0.0);

  ref_solver->IntVk_phi_unc_g = ref_solver->cdf_phi_unc_group_cell; //[g][k] copy

  auto& unc_fem_tally = uncollided_pwl_tally.tally_global;

  //============================================= Determine integrals
  auto& fv_tally_block = uncollided_fv_tally;
  for (int g=0; g<G; ++g)
  {
    for (auto& cell : grid->local_cells)
    {
      auto fv_view  = ref_solver->fv->MapFeView(cell.local_id);
      auto pwl_view =  ref_solver->pwl->GetCellMappingFE(cell.local_id);
      int k = cell.local_id;

      int mat_id = cell.material_id;
      int xs_id = ref_solver->matid_xs_map[mat_id];

      auto mat = chi_physics_handler.material_stack[mat_id];
      auto xs = std::static_pointer_cast<chi_physics::TransportCrossSections>(
        mat->properties[xs_id]);

      double siga = xs->sigma_ag[0];
      double sigt = xs->sigma_tg[0];
      double sigs = sigt-siga;

      const int num_particles = 1000;
      double sum_of_abs_point_vals = 0.0;
      std::vector<double> shape_values;

      for (int i=0; i<num_particles; ++i)
      {
        auto position = GetRandomPositionInCell(&ref_solver->rng0, cell_vol_info[k]);

        pwl_view->ShapeValues(position,shape_values);

        for (int dof=0; dof<pwl_view->num_nodes; ++dof)
        {
          int irfem = ref_solver->pwl->MapDOFLocal(cell, dof, ref_solver->uk_man_pwld,/*m*/0,/*g*/0);
          sum_of_abs_point_vals +=
            std::fabs(shape_values[dof] * unc_fem_tally[irfem]);
        }//for dof
      }//for i
      double avg_abs_value = sigs*sum_of_abs_point_vals/num_particles;
      double IntVk_phi_g_val = sigs* avg_abs_value * fv_view->volume;

      ref_solver->IntVSumG_phi_unc      += IntVk_phi_g_val;
      ref_solver->IntV_phi_unc_g[g]     += IntVk_phi_g_val;
      ref_solver->IntVk_phi_unc_g[g][k] += IntVk_phi_g_val;
    }//for lc
  }//for g

  //============================================= Build group cdf
  double sum_g_IntV_phi_unc_g = 0.0;
  for (int g=0; g<G; ++g)
  {
    sum_g_IntV_phi_unc_g += ref_solver->IntV_phi_unc_g[g];
    ref_solver->cdf_phi_unc_group[g] = sum_g_IntV_phi_unc_g /
                                       ref_solver->IntVSumG_phi_unc;
  }

  //============================================= Build group cell cdf
  for (int g=0; g<G; ++g)
  {
    double sum_k_IntVk_phi_unc_g = 0.0;
    for (int k=0; k<num_local_cells; ++k)
    {
      sum_k_IntVk_phi_unc_g += ref_solver->IntVk_phi_unc_g[g][k];
      ref_solver->cdf_phi_unc_group_cell[g][k] = sum_k_IntVk_phi_unc_g /
                                                 ref_solver->IntV_phi_unc_g[g];
    }//for k
  }//for g

  chi_log.Log(LOG_ALL) << "IntVSumG_phi_unc: " << ref_solver->IntVSumG_phi_unc;
}