#ifndef CHI_MONTECARLON_H
#define CHI_MONTECARLON_H

#include <map>
#include "ChiMesh/chi_mesh.h"
#include "ChiMesh/chi_meshvector.h"

namespace chi_math
{
  class RandomNumberGenerator;
}

//######################################################### Namespace
namespace mcpartra
{
  class SourceBase;
  class BoundarySource;
  class MaterialSource;
  class ResidualSourceA;
  class ResidualSourceB;

  class Solver;

  struct Particle;

  chi_mesh::Vector3 SampleRandomDirection(chi_math::RandomNumberGenerator& rng);
  std::pair<double,double> OmegaToPhiThetaSafe(const chi_mesh::Vector3& omega);
  size_t SampleCDF(const std::vector<double>& cdf, chi_math::RandomNumberGenerator& rng);
}

#endif //CHI_MONTECARLON_H