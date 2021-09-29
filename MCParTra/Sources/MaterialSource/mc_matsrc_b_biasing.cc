#include "mc_material_source.h"

#include "SourceDrivenSolver/sdsolver.h"

void mcpartra::MaterialSource::BiasCDFs(bool apply)
{

  //======================================== Initialize biased cdfs and corrections
  group_element_biased_cdf      = group_element_cdf;
  group_element_biased_cdf_corr = group_element_cdf;
  group_biased_cdf              = group_cdf;

  for (auto& element_cdf_corr : group_element_biased_cdf_corr)
    element_cdf_corr.assign(element_cdf_corr.size(), 1.0);

  if (not apply) return;

  const auto& importances = ref_solver.local_cell_importance;

  //======================================== Create biased unnormalized PDF
  double IntV_Q_total_local = 0.0;
  std::vector<double> IntV_Q_g(num_groups, 0.0);
  std::vector<GrpSrc> group_sources_biased = group_sources;
  {
    size_t g = 0;
    for (auto& group_src : group_sources_biased)
    {
      size_t elem = 0;
      for (auto& src_elem_pair : group_src)
      {
        const auto& element = src_elem_pair.second;
        const uint64_t cell_local_id = element.ParentCellLocalID();

        uint64_t dof_map = cell_local_id * num_groups + g;

        const double importance = importances[dof_map];

        src_elem_pair.first *= importance;

        IntV_Q_g[g] += src_elem_pair.first;

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
        group_sources_biased[g][elem].first / IntV_Q_g[g];
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
    for (auto& src_element : group_source)
    {
      elem_running_total += src_element.first;
      group_element_biased_cdf[g][elem] = (IntV_Q_g[g] > 0.0)?
        elem_running_total/IntV_Q_g[g] : 0.0;
      ++elem;
    }
    ++g;
  }
}