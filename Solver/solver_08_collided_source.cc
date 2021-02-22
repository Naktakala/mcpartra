#include "solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog& chi_log;

#include <chi_mpi.h>
extern ChiMPI& chi_mpi;

#include <ChiPhysics/chi_physics.h>
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Develop the collided source.*/
void chi_montecarlon::Solver::
  DevelopCollidedSource(chi_montecarlon::GridTallyBlock& input_fv_tally)
{
//  size_t num_local_cells = grid->local_cell_glob_indices.size();
//
//  //============================================= Initialize sizes and integrals
//  cdf_phi_unc_group.assign(num_grps, 0.0);
//  cdf_phi_unc_group_cell.resize(num_grps);
//  for (int g=0; g<num_grps; ++g)
//    cdf_phi_unc_group_cell[g].assign(num_local_cells, 0.0);
//  IntVSumG_phi_unc = 0.0;
//  IntV_phi_unc_g.assign(num_grps,0.0);
//
//  auto IntVk_phi_unc_g = cdf_phi_unc_group_cell; //[g][k] copy
//
//  //============================================= Determine integrals
//  auto& fv_tally_block = input_fv_tally;
//  for (int g=0; g<num_grps; ++g)
//  {
//    for (auto& cell : grid->local_cells)
//    {
//      auto fv_view = fv->MapFeView(cell.local_id);
//      int k = cell.local_id;
//
//      int mat_id = cell.material_id;
//      int xs_id = matid_xs_map[mat_id];
//
//      chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];
//      auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];
//
//      double siga = xs->sigma_ag[0];
//      double sigt = xs->sigma_tg[0];
//      double sigs = sigt-siga;
//
//      int ir = fv->MapDOFLocal(&cell, &uk_man_fv, 0/*m*/, g);
//
//      double IntVk_phi_g_val = sigs*
//        std::fabs(fv_tally_block.tally_global[ir] * fv_view->volume);
//
//      IntVSumG_phi_unc      += IntVk_phi_g_val;
//      IntV_phi_unc_g[g]     += IntVk_phi_g_val;
//      IntVk_phi_unc_g[g][k] += IntVk_phi_g_val;
//    }//for lc
//  }//for g
//
//  //============================================= Build group cdf
//  double sum_g_IntV_phi_unc_g = 0.0;
//  for (int g=0; g<num_grps; ++g)
//  {
//    sum_g_IntV_phi_unc_g += IntV_phi_unc_g[g];
//    cdf_phi_unc_group[g] = sum_g_IntV_phi_unc_g / IntVSumG_phi_unc;
//  }
//
//  //============================================= Build group cell cdf
//  for (int g=0; g<num_grps; ++g)
//  {
//    double sum_k_IntVk_phi_unc_g = 0.0;
//    for (int k=0; k<num_local_cells; ++k)
//    {
//      sum_k_IntVk_phi_unc_g += IntVk_phi_unc_g[g][k];
//      cdf_phi_unc_group_cell[g][k] = sum_k_IntVk_phi_unc_g / IntV_phi_unc_g[g];
//    }//for k
//  }//for g
//
//  chi_log.Log(LOG_ALL) << "IntVSumG_phi_unc: " << IntVSumG_phi_unc;
}