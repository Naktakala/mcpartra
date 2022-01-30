#include "sdsolver.h"

typedef unsigned long long TULL;

//###################################################################
/**Computes the relative std dev for all the tallies.*/
void mcpartra::SourceDrivenSolver::ComputeUncertainty()
{

  /**Lambda that computes the std_dev of the average.*/
  auto ComputeGridTallyStdDev = [](MultigroupTally& grid_tally,
                                   const int64_t dof_map,
                                   const uint64_t N_c)
  {
    const double x_avg  = grid_tally.tally_global[dof_map] / N_c;
    const double x2_avg = grid_tally.tally_sqr_global[dof_map] / N_c;

    const double stddev = sqrt(std::fabs(x2_avg - x_avg*x_avg) / N_c);
    const double relative_stddev = (std::isinf(stddev/x_avg))? 0.0 :
                                   (std::isnan(stddev/x_avg))? 0.0 :
                                   stddev/x_avg;

    grid_tally.tally_sigma[dof_map]          = stddev * N_c;
    grid_tally.tally_relative_sigma[dof_map] = relative_stddev;

    return std::pair<double,double>(x_avg, stddev);
  };

  if (nps_global == 0) nps_global = 1;

  //============================================= FV Tallies
  for (int t : fv_tallies)
  {
    auto& grid_tally = grid_tally_blocks[t];
    if (not grid_tally.empty())
    {
      for (auto& cell : grid->local_cells)
      {
        for (size_t m=0; m < num_moments; ++m)
        {
          for (size_t g=0; g < num_groups; ++g)
          {
            const int64_t dof_map = fv->MapDOFLocal(cell, 0, uk_man_fv, m, g);

            ComputeGridTallyStdDev(grid_tally, dof_map, nps_global);
          }//for g
        }//for m
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= PWL Tallies
  for (int t : pwl_tallies)
  {
    auto& grid_tally = grid_tally_blocks[t];
    if (not grid_tally.empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto& cell_fe_view = pwl->GetUnitIntegrals(cell);

        for (size_t v=0; v<cell.vertex_ids.size(); ++v)
        {
          for (size_t m=0; m < num_moments; ++m)
          {
            for (size_t g=0; g < num_groups; ++g)
            {
              const int64_t dof_map = pwl->MapDOFLocal(cell, v, uk_man_pwld, m, g);

              ComputeGridTallyStdDev(grid_tally, dof_map, nps_global);
            }//for g
          }//for m
        }//for node
      }//for local cell
    }//if tally active
  }//for tallies

  //============================================= Custom tallies
  for (auto& custom_tally : custom_tallies)
  {
    auto& grid_tally = custom_tally.grid_tally;

    CustomVolumeTally::TallyFluctuationChart new_chart;
    auto chart_size = grid_tally.tally_local.size();

    new_chart.average.reserve(chart_size);
    new_chart.sigma  .reserve(chart_size);

    for (size_t m=0; m < num_moments; ++m)
    {
      for (size_t g=0; g < num_groups; ++g)
      {
        const auto dof_map = uk_man_fv.MapUnknown(m, g);

        auto x_avg_stddev = ComputeGridTallyStdDev(grid_tally, dof_map, nps_global);

        //============================= Update fluctuation chart
        const double normalization = source_normalization *
                                     options.tally_multipl_factor /
                                     custom_tally.tally_volume;

        const double normalized_avg = x_avg_stddev.first  * normalization;
        const double normalized_std = x_avg_stddev.second * normalization;

        new_chart.average.push_back(normalized_avg);
        new_chart.sigma  .push_back(normalized_std);

      }//for g
    }//for m

    custom_tally.tally_fluctation_chart.push_back(new_chart);
  }
}