#ifndef _mc_rmc2_source_h
#define _mc_rmc2_source_h

#include "../mc_base_source.h"

#include <ChiMath/Quadratures/quadrature_gausslegendre.h>

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

#include "SourceDrivenSolver/sdsolver.h"

//###################################################################
/**Residual source class.*/
class mcpartra::ResidualSourceB : public mcpartra::SourceBase
{
public:
  std::shared_ptr<chi_physics::FieldFunction> resid_ff;
public:
  int ref_bndry=-1;
  double ref_bndry_val = 0.0;

  bool ray_trace_phase = true;

  enum class CollidedSrcMode
  {
    STAGGERED = 1,
    DIRECT    = 2
  };

  CollidedSrcMode collided_source_mode = CollidedSrcMode::DIRECT;

private:
//  chi_math::QuadratureGaussLegendre quadrature;

  //cell_g_index, face_num, RotationMatrix, Area
  typedef std::tuple<int,int,chi_mesh::Matrix3x3,double> SourcePatch;
  std::vector<SourcePatch> source_patches;

  std::vector<double>      source_patch_cdf;
  double total_patch_area = 0.0;

  struct CellSideData
  {
    double volume;
    chi_mesh::Vector3 ref_point;
    std::vector<chi_mesh::Vector3> legs;
  };
  typedef std::pair<double,std::vector<CellSideData>> CellSideInfo;
  std::vector<CellSideInfo> cell_vol_info;

  const bool sample_uniformly;

  mcpartra::MultigroupTally uncollided_fv_tally;
  mcpartra::MultigroupTally uncollided_pwl_tally;

public:
  ResidualSourceB(mcpartra::SourceDrivenSolver& solver,
                  std::shared_ptr<chi_physics::FieldFunction>& in_resid_ff,
                  bool use_uniform_sampling=false,
                  double in_bndry_val=0.0) :
    SourceBase(SourceType::RESIDUAL_TYPE_B, solver),
    resid_ff(in_resid_ff),
    sample_uniformly(use_uniform_sampling),
    ref_bndry_val(in_bndry_val)
  {}

  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>&   ref_fv_sdm) override;

  void BuildCellVolInfo(chi_mesh::MeshContinuumPtr  ref_grid,
                        std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm);
  static chi_mesh::Vector3 GetRandomPositionInCell(
    chi_math::RandomNumberGenerator* rng,
          CellSideInfo& cell_side_info);

  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng)
  {
    if (ray_trace_phase)
      return CreateBndryParticle(rng);
    else
      return CreateCollidedParticle(rng);
  }

  mcpartra::Particle
  CreateBndryParticle(chi_math::RandomNumberGenerator& rng);

  mcpartra::Particle
  CreateCollidedParticle(chi_math::RandomNumberGenerator& rng);

//  double GetRMCParallelRelativeSourceWeight();

  bool CheckForReExecution() override;

  void DevelopRMCCollidedSource();
};



#endif