#ifndef _montecarlon_source_h
#define _montecarlon_source_h

#include"../chi_montecarlon.h"
#include "mc_base_source.h"
#include"../chi_montecarlon_particle.h"

#include <ChiMesh/chi_mesh.h>
#include <FiniteVolume/fv.h>

#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include <iomanip>


namespace chi_montecarlon
{
  enum SourceTypes
  {
    BASE_SRC     = 0,
    POINT_SRC    = 1,
    BNDRY_SRC    = 2,
    MATERIAL_SRC = 3,
    RESID_SRC    = 4,
    RESID_SRC_SU = 5,
    RESID_MOC    = 6,
    RESID_MOC_SU = 7,
    RESIDUAL     = 8,
    RESIDUAL_SU  = 9,
    RESIDUAL3    = 10
  };
}

//######################################################### Class Def
/**Parent Monte carlo source.
This source is an isotropic source at [0 0 0] with energy of 4 MeV.*/
class chi_montecarlon::Source
{
public:
  int particles_C;
  int particles_L;
  int particles_R;
  double weights_L;
  double weights_R;

public:
  chi_mesh::MeshContinuum* grid;
  SpatialDiscretization_FV*   fv_sdm;
  SourceTypes type_index;

  chi_montecarlon::Solver* ref_solver;

public:
  //00
          Source();
  //01
  virtual void Initialize(chi_mesh::MeshContinuum* ref_grid,
                          SpatialDiscretization_FV*   ref_fv_sdm,
                          chi_montecarlon::Solver* ref_solver);
  virtual chi_montecarlon::Particle
          CreateParticle(chi_math::RandomNumberGenerator* rng);

  virtual double GetParallelRelativeSourceWeight() {return 1.0;}

};

#endif