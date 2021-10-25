#include "mc_rmcA_source.h"

#include "SourceDrivenSolver/sdsolver.h"

void mcpartra::ResidualSourceA::BiasCDFs(bool apply)
{
  //======================================== Initialize biased cdfs and
  //                                         corrections
  group_element_biased_cdf      = group_element_cdf;
  group_element_biased_cdf_corr = group_element_cdf;
  group_biased_cdf              = group_cdf;

  for (auto& element_cdf_corr : group_element_biased_cdf_corr)
    element_cdf_corr.assign(element_cdf_corr.size(), 1.0);

  if (not apply) return;

  const auto& importance_info = ref_solver.local_cell_importance_info;

  //======================================== Make simple copy of group-source
  //                                         strenghts
  std::vector<std::vector<double>> group_sources_biased(num_groups);

  {
    size_t g = 0;
    for (auto& group_src_biased : group_sources_biased)
    {
      size_t num_elements = group_sources[g].size();
      group_src_biased.resize(num_elements, 0.0);
      size_t elem=0;
      for (const auto& element : group_sources[g])
      {
        group_src_biased[elem] = group_sources[g][elem]->Rstar_absolute;
        ++elem;
      }
      ++g;
    }//for g
  }

  //======================================== Create biased unnormalized PDF
  const size_t num_local_cells = grid->local_cells.size();
  const size_t num_cell_R_vals = num_local_cells*num_groups;
  double IntV_Q_total_local = 0.0;
  std::vector<double> IntV_Q_g(num_groups, 0.0);
  std::vector<double> R_abs_cellk_interior_biased(num_cell_R_vals, 0.0);

  {
    size_t g = 0;
    for (auto& group_src : group_sources)
    {
      size_t elem = 0;
      for (auto& residual_info : group_src)
      {
        if (residual_info->type == RessidualInfoType::Face) continue;
        const uint64_t cell_local_id = residual_info->cell_local_id;

        const uint64_t mg_info_map = cell_local_id * num_groups + g;
        const auto& cell_imp_info = importance_info[mg_info_map];

        const double importance = cell_imp_info.importance;

        group_sources_biased[g][elem] *= importance;

        IntV_Q_g[g] += group_sources_biased[g][elem];

        R_abs_cellk_interior_biased[mg_info_map] +=
          group_sources_biased[g][elem];

        ++elem;
      }//for src_elem_pair
      IntV_Q_total_local += IntV_Q_g[g];
      ++g;
    }//for group g
  }

  //======================================== Compute normalized pdf
  std::vector<std::vector<double>> group_element_biased_pdf(num_groups);
  for (size_t g=0; g<num_groups; ++g)
  {
    size_t num_elems = group_sources_biased[g].size();
    group_element_biased_pdf[g].resize(num_elems, 0.0);
    for (size_t elem=0; elem<num_elems; ++elem)
      group_element_biased_pdf[g][elem] =
        group_sources_biased[g][elem] / IntV_Q_g[g];
  }

  //======================================== Compute bias corrections
  for (size_t g=0; g<num_groups; ++g)
  {
    size_t num_elems = group_sources_biased[g].size();
    for (size_t elem=0; elem<num_elems; ++elem)
      group_element_biased_cdf_corr[g][elem] =
        group_element_pdf[g][elem] / group_element_biased_pdf[g][elem];
  }

  //======================================== Compute biased CDFs
  // Group sampling CDF
  group_biased_cdf.assign(num_groups, 0.0);
  double running_total = 0.0;
  for (size_t g=0; g<num_groups; ++g)
  {
    running_total += IntV_Q_g[g];
    group_biased_cdf[g] = (IntV_Q_total_local > 0.0)?
      running_total / IntV_Q_total_local : 0.0;
  }

  // Element sampling CDF
  size_t g = 0;
  for (const auto& group_source : group_sources_biased)
  {
    double elem_running_total = 0.0;
    size_t elem = 0;
    for (const auto& src_element : group_source)
    {
      elem_running_total += src_element;
      group_element_biased_cdf[g][elem] = (IntV_Q_g[g] > 0.0)?
        elem_running_total/IntV_Q_g[g] : 0.0;
      ++elem;
    }
    ++g;
  }

  //============================================= Export interior source
  //                                              as FieldFunction
  auto fv_sd = std::dynamic_pointer_cast<SpatialDiscretization>(fv_sdm);
  auto R_ff = std::make_shared<chi_physics::FieldFunction>(
    "R_interior_biased",                          //Text name
    fv_sd,                                        //Spatial Discretization
    &R_abs_cellk_interior_biased,                 //Data
    ref_solver.uk_man_fv,                         //Nodal variable structure
    0, 0);                                        //Reference variable and component

    R_ff->ExportToVTKFV("ZRoutBiased","R_interior_biased");
}