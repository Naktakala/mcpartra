#include "solver_montecarlon.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiPhysics/PhysicsMaterial/chi_physicsmaterial.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "chi_log.h"

extern ChiPhysics&  chi_physics_handler;
extern ChiLog& chi_log;

#include<iostream>
#include<math.h>
#include<typeinfo>

//###################################################################
/**Processes a scattering event from cross-sections.*/
std::pair<int,chi_mesh::Vector3>
mcpartra::Solver::
ProcessScattering(mcpartra::Particle &prtcl,
                  std::shared_ptr<chi_physics::TransportCrossSections> xs)
{
  //=================================== Sample energy
  int g      = prtcl.egrp;
  int gprime = xs->Sample_gprime(g,rng0.Rand());

  //=================================== Sample scattering cosine
  double mu = xs->SampleMu_gprime_g(g, gprime, rng0.Rand(), options.force_isotropic);

  double theta    = acos(mu);
  double varphi   = rng0.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  //=================================== Build rotation matrix

  chi_mesh::Matrix3x3 R;

  if ((prtcl.dir.z < ( 1.0-1.0e-12)) and
      (prtcl.dir.z > (-1.0+1.0e-12)))
  {
    chi_mesh::Vector3 tangent = prtcl.dir.Cross(chi_mesh::Vector3(0.0,0.0,1.0));
    chi_mesh::Vector3 binorm = prtcl.dir.Cross(tangent);

    R.SetColJVec(0,tangent/tangent.Norm());
    R.SetColJVec(1,binorm/binorm.Norm());
    R.SetColJVec(2,prtcl.dir);
  } else
  {
    R.SetDiagonalVec(1.0,1.0,1.0);
  }

  //=================================== Apply rotation matrix
  chi_mesh::Vector3 dirf = R*ref_dir;

  return {gprime,dirf/dirf.Norm()};
}


//###################################################################
/**Processes a change in importance by either playing a game of
 * Russian-Roulette or by splitting the particle.*/
void mcpartra::Solver::ProcessImportanceChange(Particle& prtcl)
{
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