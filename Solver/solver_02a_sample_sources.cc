#include "solver_montecarlon.h"

//###################################################################
/**Sample sources.*/
mcpartra::Particle mcpartra::Solver::
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

  return source->CreateParticle(rng);
}