#ifndef _montecarlon_bndry_source_h
#define _montecarlon_bndry_source_h

#include "../mc_base_source.h"

#include <ChiMesh/chi_meshmatrix3x3.h>
#include <ChiMath/chi_math.h>


//###################################################################
/**Boundary source class.*/
class mcpartra::BoundarySource : public mcpartra::SourceBase
{
public:
  const int ref_bndry;
private:
  //cell_g_index, face_num, RotationMatrix, Area
  typedef std::tuple<int,int,chi_mesh::Matrix3x3,double> SourcePatch;
  std::vector<SourcePatch> source_patches;

  std::vector<double>      source_patch_cdf;
public:
  BoundarySource(mcpartra::SourceDrivenSolver& solver, const int in_ref_bndry) :
    SourceBase(SourceType::BNDRY_SRC,solver),
    ref_bndry(in_ref_bndry)
  {}

  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm) override;

  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

//  double GetParallelRelativeSourceWeight() override;
};

#endif