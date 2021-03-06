#ifndef MCPARTRA_SOURCE_BASE_H
#define MCPARTRA_SOURCE_BASE_H

#include"../chi_montecarlon.h"
#include "mc_base_source.h"
#include"../chi_montecarlon_particle.h"

#include <ChiMesh/chi_mesh.h>
#include <FiniteVolume/fv.h>

#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include <iomanip>

#include <stdexcept>


namespace mcpartra
{
enum SourceType
{
  BASE_SRC = 0,
  POINT_SRC = 1,
  BNDRY_SRC = 2,
  MATERIAL_SRC = 3,
  RESIDUAL_TYPE_A = 4,
  RESIDUAL_TYPE_B = 5,
};

//######################################################### Class Def
/**Parent Monte carlo source.
This source is an isotropic source at [0 0 0] with energy of 4 MeV.*/
class SourceBase
{
private:
  SourceType type_index;
public:
  chi_mesh::MeshContinuumPtr grid = nullptr;
  std::shared_ptr<SpatialDiscretization_FV> fv_sdm = nullptr;

  mcpartra::Solver& ref_solver;

private:
  double local_source_rate = 0.0;
  double globl_source_rate = 0.0;
  bool source_rate_determined = false;

public:
  //00
  /**Constructor*/
  SourceBase(SourceType type_spec,
             mcpartra::Solver& solver) :
         type_index(type_spec),
         ref_solver(solver)
  {};

  //01
  virtual void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                          std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm);

  virtual mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng);

  virtual bool CheckForReExecution() { return false; }

  SourceType Type() const {return type_index;}

protected:
  void SetSourceRates(double in_local_source_rate,
                      double in_globl_source_rate)
  {
    local_source_rate = in_local_source_rate;
    globl_source_rate = in_globl_source_rate;
    source_rate_determined = true;
  }

public:
  double LocalSourceRate() const
  {
    if (not source_rate_determined)
      throw std::logic_error("Source rate requested without being determined.");
    return local_source_rate;
  }
  double GlobalSourceRate() const
  {
    if (not source_rate_determined)
      throw std::logic_error("Source rate requested without being determined.");
    return globl_source_rate;
  }

  virtual ~SourceBase() = default;
};
}//namespace montecarlon

#endif //MCPARTRA_SOURCE_BASE_H