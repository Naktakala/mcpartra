#include "../mcpartra.h"

#include "ChiMath/RandomNumberGeneration/random_number_generator.h"
#include "Sources/mc_volume_src_element.h"
#include "Sources/mc_surface_src_element.h"

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
/**Samples and exponential importance representation.*/
std::pair<chi_mesh::Vector3,double> mcpartra::
  SampleSpecialRandomDirection(chi_math::RandomNumberGenerator &rng,
                                 const chi_mesh::Vector3& omega_J,
                                 const std::pair<double,double>& a_b_pair)
{
  const double a = a_b_pair.first;
  const double b = a_b_pair.second;

  if (std::fabs(b) < 1.0e-8)
    return std::make_pair(SampleRandomDirection(rng),1.0);

  //======================================== Define utilities
  const double TWO_PI = 2.0*M_PI;

  auto Isotropic_PDF = [](double mu) {return 0.5;};
  auto Exponential_PDF = [TWO_PI,a,b](double mu) {return TWO_PI * exp(a + b * mu );};

  //======================================== Rejection sample pdf for mu
//  // Find domain size
//  double max_psi = 0.0;
//  max_psi = std::max(max_psi, PDF(-1.0));
//  max_psi = std::max(max_psi, PDF( 1.0));
//
//  bool rejected = true;
//  double mu_prime = 1.0;
//  for (int i=0; i<10000; ++i)
//  {
//    mu_prime = rng.Rand() * 2.0 - 1.0;
//    const double random_PDF = rng.Rand() * max_psi;
//
//    if (random_PDF < PDF(mu_prime)) rejected = false;
//    if (not rejected) break;
//  }

  const double C_0 = (TWO_PI/b) * exp(a - b);
  const double theta_dvi_C0 = rng.Rand()/C_0;
  double mu_prime = (1.0/b) * log( exp(-b) * (theta_dvi_C0 + 1) );

  double weight_correction = Isotropic_PDF(mu_prime) / Exponential_PDF(mu_prime);

  //======================================== Compute omega in ref-coordinates
  //Sample direction
  double theta  = acos(mu_prime);
  double varphi = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 omega_prime;
  omega_prime.x = sin(theta) * cos(varphi);
  omega_prime.y = sin(theta) * sin(varphi);
  omega_prime.z = cos(theta);

  //======================================== Perform rotation
  //Build rotation matrix
  chi_mesh::Matrix3x3 R;

  const chi_mesh::Vector3 khat(0.0,0.0,1.0);

  if      (omega_J.Dot(khat) >  0.9999999)
    R.SetDiagonalVec(1.0, 1.0, 1.0);
  else if (omega_J.Dot(khat) < -0.9999999)
    R.SetDiagonalVec(1.0, 1.0,-1.0);
  else
  {
    chi_mesh::Vector3 binorm = khat.Cross(omega_J);
    binorm = binorm/binorm.Norm();

    chi_mesh::Vector3 tangent = binorm.Cross(omega_J);
    tangent = tangent/tangent.Norm();

    R.SetColJVec(0, tangent);
    R.SetColJVec(1, binorm);
    R.SetColJVec(2, omega_J);
  }

  chi_mesh::Vector3 omega = R * omega_prime;

  return std::make_pair(omega, weight_correction);
}

//###################################################################
/**Gets a cosine law random direction relative to a normal.*/
chi_mesh::Vector3 mcpartra::
  RandomCosineLawDirection(chi_math::RandomNumberGenerator& rng,
                           const chi_mesh::Vector3& normal)
{
  //Build rotation matrix
  chi_mesh::Matrix3x3 R;

  chi_mesh::Vector3 khat(0.0,0.0,1.0);

  if      (normal.Dot(khat) >  0.9999999)
    R.SetDiagonalVec(1.0,1.0,1.0);
  else if (normal.Dot(khat) < -0.9999999)
    R.SetDiagonalVec(1.0,1.0,-1.0);
  else
  {
    chi_mesh::Vector3 binorm = khat.Cross(normal);
    binorm = binorm/binorm.Norm();

    chi_mesh::Vector3 tangent = binorm.Cross(normal);
    tangent = tangent/tangent.Norm();

    R.SetColJVec(0,tangent);
    R.SetColJVec(1,binorm);
    R.SetColJVec(2,normal);
  }

  //Sample direction
  double costheta = rng.Rand();     //Sample half-range only
  double theta    = acos(sqrt(costheta));
  double varphi   = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  return R*ref_dir;
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

//###################################################################
/**Samples the cell interior*/
chi_mesh::Vector3 mcpartra::
  GetRandomPositionInCell(chi_math::RandomNumberGenerator& rng,
                          const mcpartra::CellGeometryData &cell_info)
{
  chi_mesh::Vector3 position;

  const auto& cdf = cell_info.volume_elements_cdf;
  int64_t element_id = std::lower_bound(cdf.begin(),
                                        cdf.end(),
                                        rng.Rand()) - cdf.begin();

  const auto& element = cell_info.volume_elements[element_id];
  position = element.SampleRandomPosition(rng);


  return position;
}

  //###################################################################
  /**Samples the cell interior*/
  chi_mesh::Vector3 mcpartra::
  GetRandomPositionOnCellFace(
    chi_math::RandomNumberGenerator& rng,
    const mcpartra::CellGeometryData &cell_info,
    const size_t face_index,
    int* face_sampled/*=nullptr*/,
    bool random_face/*=false*/)
{
  auto face_id = static_cast<int64_t>(face_index);

  //=================================== If random face
  if (random_face)
  {
    const auto& face_cdf = cell_info.faces_cdf;
    face_id = std::lower_bound(face_cdf.begin(),
                               face_cdf.end(),
                               rng.Rand()) - face_cdf.begin();
  }

  //=================================== Specific face identified
  if (face_sampled != nullptr)
    *face_sampled = static_cast<int>(face_id);

  const auto& element_cdf = cell_info.faces_surface_elements_cdf[face_id];
  const auto element_id = std::lower_bound(element_cdf.begin(),
                                           element_cdf.end(),
                                           rng.Rand()) - element_cdf.begin();

  const auto& element = cell_info.faces_surface_elements[face_id][element_id];

  return element.SampleRandomPosition(rng);
}


//###################################################################
/**Rotates a vector v about an axis k by angle defined according to
 * the right-hand rule.*/
chi_mesh::Vector3 mcpartra::
  RotateVec3AboutAxisRHL(const chi_mesh::Vector3 &v,
                         const chi_mesh::Vector3 &k,
                         const double angle)
{
  return v*cos(angle) +
           k.Cross(v)*sin(angle) +
           k*(k.Dot(v))*(1.0-cos(angle));
}