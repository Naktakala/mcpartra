#include "../mcpartra.h"

#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include "ChiMesh/chi_meshvector.h"

//###################################################################
/**Samples a random direction and returns the direction vector.*/
chi_mesh::Vector3 mcpartra::
  SampleRandomDirection(chi_math::RandomNumberGenerator &rng)
{
  double costheta = 2.0*rng.Rand() - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng.Rand()*2.0*M_PI;

  return chi_mesh::Vector3{sin(theta) * cos(varphi),
                           sin(theta) * sin(varphi),
                           cos(theta)};
}

//###################################################################
/**Determines the azimuthal- and polar-angle associated with
 * the given direction vector.
 * Returns a pair = [azimuthal-angle,polar-angle].*/
std::pair<double,double> mcpartra::
  OmegaToPhiThetaSafe(const chi_mesh::Vector3 &omega)
{
  // Notes: asin maps [-1,+1] to [-pi/2,+pi/2]
  //        acos maps [-1,+1] to [0,pi]
  // This mapping requires some logic for determining the azimuthal angle.
  //
  const auto omega_hat = omega.Normalized();

  double mu = omega_hat.z;
  mu = std::min(mu,  1.0);
  mu = std::max(mu, -1.0);

  double theta = acos(mu);

  //=================================== Handling omega aligned to k_hat
  if (std::fabs(omega_hat.z) < 1.0e-16) return {0.0,theta};

  //=================================== Computing varphi for NE and NW quadrant
  if (omega_hat.y >= 0.0)
  {
    double cos_phi = omega_hat.y/sin(theta);
    cos_phi = std::min(cos_phi,  1.0);
    cos_phi = std::max(cos_phi, -1.0);

    double phi = acos(cos_phi);

    return {phi, theta};
  }
  //=================================== Computing varphi for SE and SW quadrant
  else
  {
    double cos_phi = omega_hat.y/sin(theta);
    cos_phi = std::min(cos_phi,  1.0);
    cos_phi = std::max(cos_phi, -1.0);

    double phi = 2.0*M_PI - acos(cos_phi);

    return {phi, theta};
  }
}

//###################################################################
/**Samples a CDF.*/
size_t mcpartra::
  SampleCDF(const std::vector<double> &cdf, chi_math::RandomNumberGenerator &rng)
{
  auto bin_iterator = std::lower_bound(cdf.begin(), cdf.end(), rng.Rand());

  return static_cast<size_t>(std::distance(cdf.begin(), bin_iterator));
}