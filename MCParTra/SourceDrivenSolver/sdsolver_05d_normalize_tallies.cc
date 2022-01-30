#include "sdsolver.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Compute tally square contributions.*/
void mcpartra::SourceDrivenSolver::NormalizeTallies()
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Normalizing tallies.";

  if (nps_global == 0) nps_global = 1;

  const double normalization = source_normalization *
                               options.tally_multipl_factor /
                               static_cast<double>(nps_global);

  //======================================== FV Tallies
  for (unsigned int t : fv_tallies)
  {
    auto& grid_tally = grid_tally_blocks[t];
    if (not grid_tally.empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto cell_fv_view = fv->MapFeView(cell.local_id);
        const double cell_volume = cell_fv_view->volume;

        for (size_t m=0; m < num_moments; ++m)
        {
          for (size_t g=0; g < num_groups; ++g)
          {
            const int64_t dof_map = fv->MapDOFLocal(cell, 0, uk_man_fv, m, g);

            grid_tally.tally_global[dof_map] *= normalization / cell_volume;
            grid_tally.tally_sigma[dof_map]  *= normalization / cell_volume;
          }//for g
        }//for m
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= PWL Tallies
  for (unsigned int t : pwl_tallies)
  {
    auto& grid_tally = grid_tally_blocks[t];
    if (not grid_tally.empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto& cell_pwl_view   = pwl->GetUnitIntegrals(cell);

        for (size_t i=0; i < cell.vertex_ids.size(); ++i)
        {
          const double shape_volume = cell_pwl_view.IntV_shapeI(i);

          for (size_t m=0; m < num_moments; ++m)
          {
            for (size_t g=0; g < num_groups; ++g)
            {
              const int64_t dof_map = pwl->MapDOFLocal(cell, i, uk_man_pwld, m, g);

              grid_tally.tally_global[dof_map] *= normalization / shape_volume;
              grid_tally.tally_sigma[dof_map]  *= normalization / shape_volume;
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
    const double tally_volume = custom_tally.tally_volume;

    for (size_t m=0; m < num_moments; ++m)
    {
      for (size_t g=0; g < num_groups; ++g)
      {
        const auto ir = uk_man_fv.MapUnknown(m, g);

        grid_tally.tally_global[ir] *= normalization / tally_volume;
        grid_tally.tally_sigma[ir]  *= normalization / tally_volume;
      }//for g
    }//for m
  }//for custom_tally

  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done normalizing tallies.";

}