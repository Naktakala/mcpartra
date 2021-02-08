#include "solver_montecarlon.h"

//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::NormalizeTallies()
{
  if (nps_global == 0) nps_global = 1;

  //======================================== FV Tallies
  for (int t : fv_tallies)
  {
    if (not grid_tally_blocks[t].empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto cell_fv_view = fv->MapFeView(cell.local_id);

        for (int m=0; m<num_moms; ++m)
        {
          for (int g=0; g<num_grps; ++g)
          {
            int ir = fv->MapDOFLocal(&cell, &dof_structure_fv, m, g);

            grid_tally_blocks[t].tally_global[ir] *=
              source_normalization *
              tally_multipl_factor/nps_global/cell_fv_view->volume;
            grid_tally_blocks[t].tally_sigma[ir] *=
              source_normalization *
              tally_multipl_factor/cell_fv_view->volume;

          }//for g
        }//for m
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= PWL Tallies
  for (int t : pwl_tallies)
  {
    if (not grid_tally_blocks[t].empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto cell_pwl_view   = pwl->MapFeViewL(cell.local_id);

        for (int i=0; i<cell.vertex_ids.size(); ++i)
        {
          for (int m=0; m<num_moms; ++m)
          {
            for (int g=0; g<num_grps; ++g)
            {
              int ir = pwl->MapDFEMDOFLocal(&cell, i, &dof_structure_fem, m, g);

              grid_tally_blocks[t].tally_global[ir] *=
                source_normalization *
                tally_multipl_factor/nps_global/cell_pwl_view->IntV_shapeI[i];

            }//for g
          }//for m
        }//for node
      }//for local cell
    }//if tally active
  }//for tallies
}