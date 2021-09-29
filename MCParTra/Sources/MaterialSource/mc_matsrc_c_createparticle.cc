#include "mc_material_source.h"

#include "SourceDrivenSolver/sdsolver.h"
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

//###################################################################
/**Creates a source particle.*/
mcpartra::Particle mcpartra::MaterialSource::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  const std::string fname{__FUNCTION__};

  Particle new_particle;

  //======================================== Sample group
  size_t g = SampleCDF(group_biased_cdf, rng);

  new_particle.egrp = static_cast<int>(g);

  if (group_element_cdf[g].empty())
  {new_particle.alive = false; return new_particle;}

  //======================================== Sample element
  size_t elem = SampleCDF(group_element_biased_cdf[g], rng);

  auto& src_element = group_sources[g][elem].second;

  //======================================== Sample position
  new_particle.pos = src_element.SampleRandomPosition(rng);

  //======================================== Sample direction
  new_particle.dir = SampleRandomDirection(rng);

  //======================================== Determine weight
  new_particle.w = 1.0 * group_element_biased_cdf_corr[g][elem];

  new_particle.cur_cell_local_id  = src_element.ParentCellLocalID();
  new_particle.cur_cell_global_id = src_element.ParentCellGlobalID();

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::RayTraceMethod::UNCOLLIDED;

  //======================================== Determine moment indices
  const auto phi_theta = OmegaToPhiThetaSafe(new_particle.dir);
  const double& phi   = phi_theta.first;
  const double& theta = phi_theta.second;

  size_t num_moments = ref_solver.num_moments;
  for (size_t m=1; m<num_moments; ++m)
  {
    const auto& ell_em = ref_solver.m_to_ell_em_map[m];
    auto ell = static_cast<unsigned int>(ell_em.first);
    int em           = ell_em.second;

    new_particle.moment_values[m] = chi_math::Ylm(ell,em,phi,theta);
  }

  return new_particle;
}