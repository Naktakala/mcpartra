#ifndef MCPARTRA_BNDRY_SOURCE_H
#define MCPARTRA_BNDRY_SOURCE_H

#include "Sources/mc_base_source.h"
#include "Sources/mc_surface_src_element.h"

#include <ChiMesh/chi_meshmatrix3x3.h>
#include <ChiMath/chi_math.h>


//###################################################################
/**Boundary source class.*/
class mcpartra::BoundaryIsotropicSource : public mcpartra::SourceBase
{
public:
  const int ref_bndry;
  const std::vector<double> mg_isotropic_strength;
private:
  //cell_g_index, face_num, RotationMatrix, Area*strength
  typedef std::tuple<int,int,chi_mesh::Matrix3x3,double> SourcePatch;

  //Surface source element and strength
  typedef std::pair<double, SurfaceSourceElement> SourceElement;

  std::vector<std::vector<SourceElement>> group_elements;

  std::vector<std::vector<double>>        group_element_cdf;

  std::vector<double>                     group_cdf;
public:
  BoundaryIsotropicSource(mcpartra::SourceDrivenSolver& solver,
                          const int in_ref_bndry,
                          const std::vector<double> in_mg_isotropic_strength) :
    SourceBase(SourceType::BNDRY_SRC,solver),
    ref_bndry(in_ref_bndry),
    mg_isotropic_strength(in_mg_isotropic_strength)
  {}

  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
                  size_t ref_num_groups,
                  const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
                  const std::vector<CellGeometryData>& ref_cell_geometry_info) override;

  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;
};

#endif //MCPARTRA_BNDRY_SOURCE_H