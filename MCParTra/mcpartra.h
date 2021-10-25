#ifndef MCPARTRA_H
#define MCPARTRA_H

#include <map>
#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/chi_meshvector.h"

namespace chi_math
{
  class RandomNumberGenerator;
  class SparseMatrix;
}

//######################################################### Namespace
namespace mcpartra
{
  class SourceBase;
  class BoundarySource;
  class MaterialSource;
  class ResidualSourceA;

  class SourceDrivenSolver;

  struct Particle;

  enum RayTraceMethod
  {
    STANDARD   = 0,
    UNCOLLIDED = 1
  };

  class VolumeElement;
  class SurfaceSourceElement;

  /**Structure to hold all of the constituents of non-curvilinear cells.*/
  struct CellGeometryData
  {
    double total_volume=0.0;
    std::vector<VolumeElement> volume_elements;
    std::vector<double> volume_elements_cdf;

    double total_area=0.0;
    std::vector<double> faces_total_area;
    std::vector<double> faces_cdf;

    std::vector<std::vector<SurfaceSourceElement>> faces_surface_elements;
    std::vector<std::vector<double>> faces_surface_elements_cdf;
  };

  struct CellImportanceInfo
  {
    double importance = 1.0;
    chi_mesh::Vector3 omega_J;
    double a = std::log(0.5);
    double b = 0.0;

    double ExpRep(const chi_mesh::Vector3& omega) const
    {
      double mu_J = omega.Dot(omega_J);
      return exp( a + b * mu_J);
    }
  };

  //mc_utils_01
  chi_mesh::Vector3 SampleRandomDirection(chi_math::RandomNumberGenerator& rng);
  std::pair<chi_mesh::Vector3,double>
    SampleSpecialRandomDirection(chi_math::RandomNumberGenerator& rng,
                                 const chi_mesh::Vector3& omega_J,
                                 const std::pair<double,double>& a_b_pair);
  chi_mesh::Vector3 RandomCosineLawDirection(chi_math::RandomNumberGenerator& rng,
                                             const chi_mesh::Vector3& normal);
  std::pair<double,double> OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega);
  size_t SampleCDF(const std::vector<double>& cdf, chi_math::RandomNumberGenerator& rng);

  chi_mesh::Vector3 GetRandomPositionInCell(chi_math::RandomNumberGenerator& rng,
                                            const CellGeometryData& cell_info);

  chi_mesh::Vector3 GetRandomPositionOnCellFace(
    chi_math::RandomNumberGenerator& rng,
    const CellGeometryData& cell_info,
    size_t face_index,
    int* face_sampled = nullptr,
    bool random_face = false);

  //mc_utils_02
  typedef std::vector<std::vector<double>> MatDbl;

  MatDbl ComputeGroupToGroupScatteringCDFs(
      size_t G,
      const chi_math::SparseMatrix& isotropic_transfer_matrix);

  typedef std::pair<double,double> MuCP; //Cosine and cumulative probability pairs
  typedef std::vector<MuCP> CosineCDF;   //CDF of cosines
  typedef std::vector<CosineCDF> VecCosineCDFs; //Vector of CDFs

  std::vector<VecCosineCDFs> ComputeDiscreteScatteringAngles(
    size_t G,
    const std::vector<chi_math::SparseMatrix>& transfer_matrices,
    size_t num_moments_to_support);
}

#endif //MCPARTRA_H