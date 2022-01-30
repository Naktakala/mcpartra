#include "mg_scattering_cdfs.h"

#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "mcpartra.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include <algorithm>
#include <cmath>
#include <cassert>

mcpartra::MultigroupScatteringCDFs::
  MultigroupScatteringCDFs(const chi_physics::TransportCrossSections& mgxs,
                           unsigned int L_parameter)
{
  chi_log.Log() << "Building mg scattering pdfs. L=" << L_parameter; //TODO: Remove
  chi_log.Log() << "  S-matrices available: " << mgxs.transfer_matrices.size();

  size_t num_groups = mgxs.num_groups;

  assert(not mgxs.transfer_matrices.empty());
  const auto& S0 = mgxs.transfer_matrices[0];

  cdf_gprime_g =
    mcpartra::ComputeGroupToGroupScatteringCDFs(num_groups, S0);

  scat_angles_gprime_g =
    mcpartra::ComputeDiscreteScatteringAngles(num_groups,
                                              mgxs.transfer_matrices,
                                              L_parameter);
}

unsigned int mcpartra::MultigroupScatteringCDFs::
  Sample_gprime(unsigned int gp, double rn) const
{
  unsigned int gto =
    std::lower_bound(cdf_gprime_g[gp].begin(),
                     cdf_gprime_g[gp].end(),
                     rn) - cdf_gprime_g[gp].begin();

  return gto;
}

double mcpartra::MultigroupScatteringCDFs::
  SampleMu_gprime_g(unsigned int gp,
                    unsigned int g, double rn, bool isotropic) const
{
  auto& chi_log = ChiLog::GetInstance();

  double mu;

  struct
  {
    bool operator()(const std::pair<double,double>& left, double val)
    {return left.second <= val;}
  }compare;

  if (isotropic or scat_angles_gprime_g[gp][g].empty())
    mu = 2.0*rn-1.0;
  else
  {
    unsigned int angle_num =
      std::lower_bound(scat_angles_gprime_g[gp][g].begin(),
                       scat_angles_gprime_g[gp][g].end(),rn,
                       compare) -
                       scat_angles_gprime_g[gp][g].begin();
    mu = scat_angles_gprime_g[gp][g][angle_num].first;
  }

  if (std::isnan(mu))
  {
    chi_log.Log(LOG_ALLERROR)
    << "Mu corruption in GolubFischer sample mu.\n"
    << " g      =" << g << "\n"
    << " gprime =" << gp << "\n"
    << " rn     =" << rn;
  }

  if (mu > 1.0)  mu = 1.0;
  if (mu < -1.0) mu = -1.0;

  return mu;
}