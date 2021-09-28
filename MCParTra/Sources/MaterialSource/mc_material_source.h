#ifndef MCPARTRA_MATERIAL_SOURCE_H
#define MCPARTRA_MATERIAL_SOURCE_H

#include "Sources/mc_base_source.h"
#include "Sources/mc_volume_src_element.h"

#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

//###################################################################
/**Material source class.*/
class mcpartra::MaterialSource : public mcpartra::SourceBase
{
private:
  typedef std::pair<double, VolumeSourceElement&> ElementSrc;
  typedef std::vector<ElementSrc> GrpSrc;
  typedef chi_physics::IsotropicMultiGrpSource IsoMGSrc;

  std::map<uint64_t, std::vector<VolumeSourceElement>> cell_elements;

  std::vector<bool> matid_has_q_flags;
  std::map<int, std::shared_ptr<IsoMGSrc>> matid_q_map;

private: //CDFs
  std::vector<double> group_cdf;

  std::vector<GrpSrc> group_sources;
  std::vector<std::vector<double>> group_element_cdf;

public:
  explicit
  MaterialSource(mcpartra::SourceDrivenSolver& solver) :
    SourceBase(SourceType::MATERIAL_SRC, solver)
  {};

  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
                  size_t ref_num_groups,
                  const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
                  const std::vector<CellGeometryData>& ref_cell_geometry_info) override;

  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;
};

#endif //MCPARTRA_MATERIAL_SOURCE_H