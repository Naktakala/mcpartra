#include "solver_montecarlon.h"

typedef unsigned long long TULL;

//###################################################################
/**Computes the relative std dev for all the tallies.*/
void chi_montecarlon::Solver::ComputeUncertainty()
{
  max_sigma = 0.0;
  max_relative_sigma = 0.0;

  max_fem_sigma = 0.0;
  max_fem_relative_sigma = 0.0;


  //============================================= FV Tallies
  double IntV_sigma = 0.0;
  double Vtot = 0.0;
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

            TULL divisor = (nps_global==0)? 1 : nps_global;

            double x_avg  = grid_tally_blocks[t].tally_global[ir]/divisor;
            double x2_avg = grid_tally_blocks[t].tally_sqr_global[ir]/divisor;

            double stddev = sqrt(std::fabs(x2_avg - x_avg*x_avg)/divisor);

            grid_tally_blocks[t].tally_sigma[ir] = stddev;

            max_sigma = std::max(max_sigma,stddev);

            if (!std::isinf(stddev/x_avg) and !std::isnan(stddev/x_avg))
            {
              grid_tally_blocks[t].tally_relative_sigma[ir] = stddev/x_avg;
              max_relative_sigma = std::max(max_relative_sigma,stddev/x_avg);
            }
            else
              grid_tally_blocks[t].tally_relative_sigma[ir] = 0.0;

            IntV_sigma +=
              grid_tally_blocks[t].tally_sigma[ir];

            Vtot += cell_fv_view->volume;
          }//for g
        }//for m
      }//for local cell
    }//if tally active
  }//for tallies
  avg_sigma = IntV_sigma / Vtot;

  //============================================= PWL Tallies
  double IntV_fem_sigma = 0.0;
  double Vtot_fem = 0.0;
  for (int t : pwl_tallies)
  {
    if (not grid_tally_blocks[t].empty())
    {
      for (auto& cell : grid->local_cells)
      {
        auto cell_fe_view = pwl->MapFeViewL(cell.local_id);

        for (int v=0; v<cell.vertex_ids.size(); ++v)
        {
          for (int m=0; m<num_moms; ++m)
          {
            for (int g=0; g<num_grps; ++g)
            {
              int ir = pwl->MapDFEMDOFLocal(&cell, v, &dof_structure_fem, m, g);

              TULL divisor = (nps_global==0)? 1 : nps_global;

              double x_avg  = grid_tally_blocks[t].tally_global[ir]/divisor;
              double x2_avg = grid_tally_blocks[t].tally_sqr_global[ir]/divisor;

              double stddev = sqrt(std::fabs(x2_avg - x_avg*x_avg)/divisor);

              grid_tally_blocks[t].tally_sigma[ir] = stddev;

              max_fem_sigma = std::max(max_fem_sigma,stddev);

              if (!std::isinf(stddev/x_avg) and !std::isnan(stddev/x_avg))
              {
                grid_tally_blocks[t].tally_relative_sigma[ir] = stddev/x_avg;
                max_fem_relative_sigma = std::max(max_relative_sigma,stddev/x_avg);
              }
              else
                grid_tally_blocks[t].tally_relative_sigma[ir] = 0.0;

              IntV_fem_sigma +=
                grid_tally_blocks[t].tally_sigma[ir]*
                  cell_fe_view->IntV_shapeI[v];
              Vtot_fem += cell_fe_view->IntV_shapeI[v];

            }//for g
          }//for m
        }//for node
      }//for local cell
    }//if tally active
  }//for tallies
  avg_fem_sigma = IntV_fem_sigma / Vtot_fem;

  //============================================= Custom tallies
  for (auto& custom_tally : custom_tallies)
  {
    auto& grid_tally = custom_tally.grid_tally;

    CustomVolumeTally::TallyFluctuationChart new_chart;
    auto chart_size = grid_tally.tally_local.size();

    new_chart.average.reserve(chart_size);
    new_chart.sigma  .reserve(chart_size);

    for (int m=0; m<num_moms; ++m)
    {
      for (int g=0; g<num_grps; ++g)
      {
        auto ir = dof_structure_fv.MapVariable(m,g);

        TULL divisor = (nps_global==0)? 1 : nps_global;

        double x_avg  = grid_tally.tally_global[ir]/divisor;
        double x2_avg = grid_tally.tally_sqr_global[ir]/divisor;

        double stddev = sqrt(std::fabs(x2_avg - x_avg*x_avg)/divisor);

        grid_tally.tally_sigma[ir] = stddev;

        if (!std::isinf(stddev/x_avg) and !std::isnan(stddev/x_avg))
          grid_tally.tally_relative_sigma[ir] = stddev/x_avg;
        else
          grid_tally.tally_relative_sigma[ir] = 0.0;

        double normalized_avg = x_avg*
                                source_normalization*
                                tally_multipl_factor/
                                custom_tally.tally_volume;

        double normalized_std = stddev*
                                source_normalization*
                                tally_multipl_factor/
                                custom_tally.tally_volume;

        new_chart.average.push_back(normalized_avg);
        new_chart.sigma  .push_back(normalized_std);

      }//for g
    }//for m

    custom_tally.tally_fluctation_chart.push_back(new_chart);
  }
}