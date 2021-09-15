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
  class ResidualSourceB;

  class SourceDrivenSolver;

  struct Particle;

  enum RayTraceMethod
  {
    STANDARD   = 0,
    UNCOLLIDED = 1
  };

  struct TallyMethod
  {
    static const int STANDARD     = 0;
    static const int RMC_CHAR_RAY = 1;
  };

  chi_mesh::Vector3 SampleRandomDirection(chi_math::RandomNumberGenerator& rng);
  std::pair<double,double> OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega);
  size_t SampleCDF(const std::vector<double>& cdf, chi_math::RandomNumberGenerator& rng);
  void ComputeDiscreteScatteringAngles(size_t G,
                                       const std::vector<chi_math::SparseMatrix>& transfer_matrices);
}

#endif //MCPARTRA_H