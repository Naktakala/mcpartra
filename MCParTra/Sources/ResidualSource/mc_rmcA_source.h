#ifndef MCPARTRA_RMC_SOURCE_A_H
#define MCPARTRA_RMC_SOURCE_A_H

#include "Sources/mc_base_source.h"
#include "Sources/mc_volume_src_element.h"
#include "Sources/mc_surface_src_element.h"

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

//###################################################################
/**Residual source class.*/
class mcpartra::ResidualSourceA : public mcpartra::SourceBase
{
public:
  std::shared_ptr<chi_physics::FieldFunction> resid_ff;
private:
  typedef std::pair<double, VolumeSourceElement&> ElementSrc;
  typedef std::vector<ElementSrc> GrpSrc;

  /**Simplified material structure that can be passed around.*/
  struct MaterialData
  {
    double siga=0.0;
    double Q=0.0;
  };

  std::unique_ptr<std::vector<CellGeometryData>> cell_geometry_info;

  /**Structure to store cell interior residual info.*/
  struct RCellInterior
  {
    uint64_t cell_local_id=0;
    double average_rstar_absolute=0.0;
    double maximum_rstar_absolute=-1.0e32;
    double Rstar_absolute=0.0;
  };
  typedef std::vector<RCellInterior> VecRCellInterior;
  std::vector<VecRCellInterior> residual_info_cell_interiors;

  /**Structure to store cell-face pairs.*/
  struct RCellFace
  {
    uint64_t cell_local_id=0;
    int ass_face=0;
    double maximum_rstar_absolute=-1.0e32;
    double Rstar_absolute=0.0;
  };
  typedef std::vector<RCellFace> VecRCellFace;
  std::vector<VecRCellFace> residual_info_cell_bndry_faces;

private: //CDFs
  std::vector<double> group_cdf;

  std::vector<GrpSrc> group_sources;
  std::vector<std::vector<double>> group_element_cdf;

private:
  std::vector<double> R_abs_cellk_interior; //for field function

  std::vector<double> R_abs_localdomain_interior; ///< Per group
  std::vector<double> R_abs_localdomain_surface;  ///< Per group

  std::vector<std::vector<double>> domain_cdf; ///< Per group
  std::vector<std::vector<double>> surface_cdf; ///< Per group

public:
  explicit
  ResidualSourceA(mcpartra::SourceDrivenSolver& solver,
                  std::shared_ptr<chi_physics::FieldFunction>& in_resid_ff) :
    SourceBase(SourceType::RESIDUAL_TYPE_A, solver),
    resid_ff(in_resid_ff)
  {}

  //a
  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
                  size_t ref_num_groups,
                  const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
                  const std::vector<CellGeometryData>& ref_cell_geometry_info) override;

  //b
  void PopulateMaterialData(int mat_id, int group_g,
                            MaterialData& mat_data);

  double GetResidualFFPhi(std::vector<double> &N_in,
                          size_t dofs,
                          uint64_t cell_local_id,
                          int egrp);

  chi_mesh::Vector3 GetResidualFFGradPhi(
                          std::vector<chi_mesh::Vector3>& Grad_in,
                          size_t dofs,
                          uint64_t cell_local_id,
                          int egrp);

  //c
  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

  //d
  void RemoveFFDiscontinuities();

}; //MCPARTRA_RMC_SOURCE_A_H



#endif