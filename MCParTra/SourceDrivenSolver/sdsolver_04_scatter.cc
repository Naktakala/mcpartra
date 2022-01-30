#include "sdsolver.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "chi_log.h"

extern ChiPhysics&  chi_physics_handler;
extern ChiLog& chi_log;

#include<math.h>

//###################################################################
/**Processes a scattering event from cross-sections.*/
std::pair<int,chi_mesh::Vector3>
mcpartra::SourceDrivenSolver::
ProcessScattering(mcpartra::Particle &prtcl,
                  const MultigroupScatteringCDFs& xs)
{
  typedef chi_mesh::Vector3 Vec3;
  const Vec3 khat(0.0,0.0,1.0);

  //=================================== Sample energy
  int g      = prtcl.egrp;
  int gprime = xs.Sample_gprime(g,rng0.Rand());

  //=================================== Sample scattering cosine
  double mu = xs.SampleMu_gprime_g(g, gprime, rng0.Rand(), options.force_isotropic);

  const double theta_L_prime    = acos(mu);
  const double varphi_L_prime   = rng0.Rand() * 2.0 * M_PI;

  const auto  n = prtcl.dir.Normalized();
  const auto  t = ((1.0 - std::fabs(n.Dot(khat)))<1.0e-8)?
                  Vec3(1,0,0) :
                  n.Cross(khat).Normalized();

  const auto dirpolar = RotateVec3AboutAxisRHL(n, t, -theta_L_prime);
  const auto dirfinal = RotateVec3AboutAxisRHL(dirpolar, n, varphi_L_prime);

  for (size_t m=0; m<num_moments; ++m)
  {
    auto& ell_em = m_to_ell_em_map[m];
    auto ell = ell_em.first;
    auto em  = ell_em.second;

    const auto angles = OmegaToPhiThetaSafe(dirfinal);
    const double varphi_L = angles.first;
    const double theta_L = angles.second;

    prtcl.moment_values[m] = chi_math::Ylm(ell, em, varphi_L, theta_L);
  }

  return {gprime,dirfinal};
}


//###################################################################
/**Processes a change in importance by either playing a game of
 * Russian-Roulette or by splitting the particle.*/
void mcpartra::SourceDrivenSolver::
  ProcessImportanceChange(Particle& prtcl, double current_cell_importance)
{
  prtcl.cur_cell_importance = current_cell_importance;
  prtcl.pre_cell_importance = prtcl.cur_cell_importance;
  if (prtcl.cur_cell_global_id == prtcl.pre_cell_global_id) return;


  //=================================== Particle splitting
  if (prtcl.cur_cell_importance > prtcl.pre_cell_importance)
  {
    double R = prtcl.pre_cell_importance/prtcl.cur_cell_importance;

    if (rng0.Rand() > R)
    {
      prtcl.pre_cell_importance = prtcl.cur_cell_importance; //prevent resplitting
      prtcl.w *= 0.5;
      Particle prtcl_split = prtcl;
      particle_source_bank.push_back(prtcl_split);
      particle_source_bank.push_back(prtcl);
      prtcl.banked = true;
      prtcl.alive = false;
    }
  }
  //=================================== Russian-Roulette
  else if (prtcl.cur_cell_importance < prtcl.pre_cell_importance)
  {
    double R = prtcl.cur_cell_importance/prtcl.pre_cell_importance;

    if (rng0.Rand() > R)
      prtcl.alive = false;
  }
}