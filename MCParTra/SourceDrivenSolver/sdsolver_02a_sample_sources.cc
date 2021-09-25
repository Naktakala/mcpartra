#include "sdsolver.h"

//###################################################################
/**Sample sources.*/
mcpartra::Particle mcpartra::SourceDrivenSolver::
  SampleSources(chi_math::RandomNumberGenerator &rng)
{
  size_t src_index = std::lower_bound(
                   local_source_cdf.begin(),
                   local_source_cdf.end(),
                   rng.Rand()) - local_source_cdf.begin();

  if (src_index >= local_source_cdf.size())
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Error sampling local_source_cdf.");

  auto& source = sources[src_index];

  auto prtcl = source->CreateParticle(rng);

  double importance = local_cell_importance[prtcl.cur_cell_local_id];
  accumulated_src_importances += importance;

  return prtcl;
}