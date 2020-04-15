#ifndef _mc_rmc2_source_h
#define _mc_rmc2_source_h

#include "../mc_base_source.h"

#include <ChiMath/Quadratures/quadrature_gausslegendre.h>

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

//###################################################################
/**Residual source class.*/
class chi_montecarlon::ResidualSource2 : public chi_montecarlon::Source
{
public:
  int ref_bndry=-1;
  double ref_bndry_val = 0.0;

  chi_physics::FieldFunction* resid_ff;
private:
  chi_math::QuadratureGaussLegendre quadrature;

  //cell_g_index, face_num, RotationMatrix, Area
  typedef std::tuple<int,int,chi_mesh::Matrix3x3,double> SourcePatch;
  std::vector<SourcePatch> source_patches;

  std::vector<double>      source_patch_cdf;

  struct CellSideData
  {
    double volume;
    chi_mesh::Vector3 ref_point;
    std::vector<chi_mesh::Vector3> legs;
  };
  typedef std::pair<double,std::vector<CellSideData>> CellSideInfo;
  std::vector<CellSideInfo> cell_vol_info;

  const bool sample_uniformly;
public:
  ResidualSource2(chi_physics::FieldFunction* in_resid_ff,
                 bool use_uniform_sampling=false,
                 double in_bndry_val=0.0);

  void Initialize(chi_mesh::MeshContinuum* ref_grid,
                  SpatialDiscretization_FV*   ref_fv_sdm,
                  chi_montecarlon::Solver* ref_solver);

  void BuildCellVolInfo(chi_mesh::MeshContinuum*  ref_grid,
                        SpatialDiscretization_FV* ref_fv_sdm);
  chi_mesh::Vector3 GetRandomPositionInCell(
          chi_montecarlon::RandomNumberGenerator* rng,
          CellSideInfo& cell_side_info);

  chi_montecarlon::Particle
  CreateBndryParticle(chi_montecarlon::RandomNumberGenerator* rng);

  chi_montecarlon::Particle
  CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);

  chi_montecarlon::Particle
  SampleBoundary(chi_montecarlon::RandomNumberGenerator* rng);

  chi_montecarlon::Particle
  DirectSampling(chi_montecarlon::RandomNumberGenerator* rng);

  double GetRMCParallelRelativeSourceWeight();
};



#endif