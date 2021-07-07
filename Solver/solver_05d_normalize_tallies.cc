#include "solver_montecarlon.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Compute tally square contributions.*/
void mcpartra::Solver::NormalizeTallies()
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Normalizing tallies.";

  if (nps_global == 0) nps_global = 1;

  //======================================== FV Tallies
  for (int _t : fv_tallies)
  {
    auto t = static_cast<unsigned int>(_t);
    if (not grid_tally_blocks[t].empty())
    {
      auto& tally_global = grid_tally_blocks[t].tally_global;
      auto& tally_sigma  = grid_tally_blocks[t].tally_sigma;
      for (auto& cell : grid->local_cells)
      {
        auto cell_fv_view = fv->MapFeView(cell.local_id);

        for (size_t m=0; m<num_moms; ++m)
        {
          for (size_t g=0; g<num_grps; ++g)
          {
            int64_t _dof_map = fv->MapDOFLocal(cell, 0, uk_man_fv, m, g);
            auto dof_map = static_cast<uint64_t>(_dof_map);

            grid_tally_blocks[t].tally_global[dof_map] *=
              source_normalization *
              options.tally_multipl_factor /
                (double)nps_global /
              cell_fv_view->volume;

            grid_tally_blocks[t].tally_sigma[dof_map] *=
              source_normalization *
              options.tally_multipl_factor /
              cell_fv_view->volume;

          }//for g
        }//for m
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= PWL Tallies
  for (int _t : pwl_tallies)
  {
    auto t = static_cast<unsigned int>(_t);
    if (not grid_tally_blocks[t].empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto& cell_pwl_view   = pwl->GetUnitIntegrals(cell);

        for (size_t i=0; i < cell.vertex_ids.size(); ++i)
        {
          for (size_t m=0; m < num_moms; ++m)
          {
            for (size_t g=0; g < num_grps; ++g)
            {
              int64_t _dof_map = pwl->MapDOFLocal(cell, i, uk_man_pwld, m, g);
              auto dof_map = static_cast<uint64_t>(_dof_map);

              grid_tally_blocks[t].tally_global[dof_map] *=
                source_normalization *
                options.tally_multipl_factor /
                (double)nps_global /
                cell_pwl_view.IntV_shapeI(i);

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
    for (size_t m=0; m<num_moms; ++m)
    {
      for (size_t g=0; g<num_grps; ++g)
      {
        auto ir = uk_man_fv.MapUnknown(m, g);

        grid_tally.tally_global[ir] *=
          source_normalization *
          options.tally_multipl_factor /
          (double)nps_global /
          custom_tally.tally_volume;

        grid_tally.tally_sigma[ir] *=
          source_normalization *
          options.tally_multipl_factor /
          custom_tally.tally_volume;

      }//for g
    }//for m
  }//for custom_tally

  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done normalizing tallies.";

}