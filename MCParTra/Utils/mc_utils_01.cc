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