#include "mc_base_source.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//#########################################################
/**Default constructor*/
chi_montecarlon::Source::Source()
{
  type_index = SourceTypes::BASE_SRC;
}

void chi_montecarlon::Source::
  Initialize(chi_mesh::MeshContinuumPtr ref_grid,
             std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm,
             chi_montecarlon::Solver* ref_solver)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
}

//#########################################################
/**Create a default particle from a point.*/
chi_montecarlon::Particle chi_montecarlon::Source::CreateParticle(chi_math::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  new_particle.pos.x = 0.0;
  new_particle.pos.y = 0.0;
  new_particle.pos.z = 0.0;

  double costheta = rng->Rand()*2.0 - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  new_particle.dir.x = sin(theta)*cos(varphi);
  new_particle.dir.y = sin(theta)*sin(varphi);
  new_particle.dir.z = cos(theta);

  new_particle.egrp = 0;
  new_particle.w = 1.0;

  chi_log.Log(LOG_0) << "Default particle created.";

  return new_particle;
}