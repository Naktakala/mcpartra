#include "solver_montecarlon.h"

//###################################################################
/**Computes PWLD transformations of the PWLD tallies.*/
void chi_montecarlon::Solver::ComputePWLDTransformations()
{
  size_t num_cells = grid->local_cell_glob_indices.size();
  if (make_pwld)
  {
    for (size_t lc=0; lc<num_cells; lc++)
    {
      int cell_g_index = grid->local_cell_glob_indices[lc];
      int map = local_cell_pwl_dof_array_address[lc];

      auto cell_pwl_view = pwl_discretization->MapFeView(cell_g_index);

      MatDbl A(cell_pwl_view->IntV_shapeI_shapeJ);
      MatDbl Ainv = chi_math::Inverse(A);
      VecDbl b(cell_pwl_view->dofs,0.0);

      for (int g=0; g<num_grps; g++)
      {
        for (int dof=0; dof<cell_pwl_view->dofs; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          b[dof] = phi_pwl_global[ir]*cell_pwl_view->IntV_shapeI[dof];
        }
        VecDbl x = chi_math::MatMul(Ainv,b);
        for (int dof=0; dof<cell_pwl_view->dofs; dof++)
        {
          int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
          phi_pwl_global[ir] = x[dof];
        }
      }
    }//for local cell lc
  }//if make_pwld
}