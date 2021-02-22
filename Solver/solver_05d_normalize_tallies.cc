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
            int ir = fv->MapDOFLocal(cell, uk_man_fv, m, g);

            grid_tally_blocks[t].tally_global[ir] *=
              source_normalization *
              tally_multipl_factor /
                (double)nps_global /
              cell_fv_view->volume;

            grid_tally_blocks[t].tally_sigma[ir] *=
              source_normalization *
              tally_multipl_factor /
              cell_fv_view->volume;

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
        auto& cell_pwl_view   = pwl->GetUnitIntegrals(cell);

        for (int i=0; i<cell.vertex_ids.size(); ++i)
        {
          for (int m=0; m<num_moms; ++m)
          {
            for (int g=0; g<num_grps; ++g)
            {
              int ir = pwl->MapDOFLocal(cell, i, uk_man_pwld, m, g);

              grid_tally_blocks[t].tally_global[ir] *=
                source_normalization *
                tally_multipl_factor /
                (double)nps_global /
                cell_pwl_view.IntV_shapeI[i];

            }//for g
          }//for m
        }//for node
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= Custom Tallies
  for (auto& custom_tally : custom_tallies)
  {
    auto& grid_tally = custom_tally.grid_tally;

    for (int m=0; m<num_moms; ++m)
    {
      for (int g=0; g<num_grps; ++g)
      {
        auto ir = uk_man_fv.MapUnknown(m, g);

        grid_tally.tally_global[ir] *=
          source_normalization *
          tally_multipl_factor /
          (double)nps_global /
          custom_tally.tally_volume;

        grid_tally.tally_sigma[ir] *=
          source_normalization *
          tally_multipl_factor /
          custom_tally.tally_volume;

      }//for g
    }//for m

  }


}