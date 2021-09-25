#include "mc_base_source.h"

#include "chi_log.h"
extern ChiLog& chi_log;

void mcpartra::SourceBase::
  Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;
}

//#########################################################
/**Create a default particle from a point.*/
mcpartra::Particle mcpartra::SourceBase::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  return Particle::MakeDeadParticle();
}