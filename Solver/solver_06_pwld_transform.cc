#include "solver_montecarlon.h"

//###################################################################
/**Computes PWLD transformations of the PWLD tallies.*/
void chi_montecarlon::Solver::ComputePWLDTransformations()
{
  if (make_pwld)
  {
    for (auto& cell : grid->local_cells)
    {
      auto cell_pwl_view = pwl->MapFeViewL(cell.local_id);

      MatDbl A(cell_pwl_view->IntV_shapeI_shapeJ);
      MatDbl Ainv = chi_math::Inverse(A);
      VecDbl b(cell_pwl_view->dofs,0.0);

//      for (int t : pwl_tallies)
      auto& raw_tally = grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]];
      auto& out_tally = grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]];

      {
        if (raw_tally.empty()) continue;

        for (int m=0; m<num_moms; ++m)
        {
          for (int g=0; g<num_grps; ++g)
          {

            for (int i=0; i<cell.vertex_ids.size(); ++i)
            {
              int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man_fem,/*m*/0,g);
              b[i] = raw_tally.tally_global[ir]*
                     cell_pwl_view->IntV_shapeI[i];
            }//for dof

            VecDbl x = chi_math::MatMul(Ainv,b);

            for (int i=0; i<cell.vertex_ids.size(); ++i)
            {
              int ir = pwl->MapDFEMDOFLocal(&cell,i,&uk_man_fem,/*m*/0,g);
              out_tally.tally_global[ir] = x[i];
            }//for dof

          }//for group
        }//for moment
      }//for tally

    }//for local cell lc
  }//if make_pwld
}