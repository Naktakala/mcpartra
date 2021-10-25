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
  /**Simplified material structure that can be passed around.*/
  struct MaterialData
  {
    double siga=0.0;
    double Q=0.0;
  };

  std::unique_ptr<std::vector<CellGeometryData>> cell_geometry_info;

  enum class RessidualInfoType : int
  {
    Interior = 0,
    Face     = 1
  };

  /**Structure to store cell interior residual info.*/
  struct RCellInterior
  {
    RessidualInfoType type = RessidualInfoType::Interior;
    uint64_t          cell_local_id=0;
    double            maximum_rstar_absolute=-1.0e32;
    double            maximum_rstar_psi_star_absolute=-1.0e32;
    double            Rstar_absolute=0.0;
    double            Rstar_psi_star_absolute=0.0;
  };

  /**Structure to store cell-face pairs.*/
  struct RCellFace : public RCellInterior
  {
    unsigned int ass_face=0;
    RCellFace() { type = RessidualInfoType::Face;}
  };

private://PDFs
  typedef std::unique_ptr<RCellInterior> RCellInteriorPtr;
  typedef std::vector<RCellInteriorPtr> GrpSrc;
  std::vector<GrpSrc>              group_sources;     ///< Per group then element
  std::vector<std::vector<double>> group_element_pdf; ///< Per group then element

private:
  std::vector<double> R_abs_cellk_interior; ///< For field function

private: //CDFs
  std::vector<std::vector<double>> group_element_cdf; ///< Per group then element
  std::vector<double>              group_cdf;         ///< Per group

private: //Biased CDFs
  std::vector<std::vector<double>> group_element_biased_cdf;      ///< Per group then element
  std::vector<std::vector<double>> group_element_biased_cdf_corr; ///< Per group then element

  std::vector<double>              group_biased_cdf;              ///< Per group
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
  void BiasCDFs(bool apply) override;

  //c
  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

  //d
  void RemoveFFDiscontinuities();

  //e
  void PopulateMaterialData(int mat_id, size_t group_g,
                            MaterialData& mat_data);

  std::vector<double>
    GetResidualFFPhiAtNodes(const chi_mesh::Cell& cell,
                            size_t num_nodes,
                            size_t variable_id,
                            size_t component_id);

  static double GetPhiH(const std::vector<double>& shape_values,
                        const std::vector<double>& phi,
                        size_t num_nodes);

  static chi_mesh::Vector3 GetGradPhiH(
    const std::vector<chi_mesh::Vector3>& grad_shape_values,
    const std::vector<double>& phi,
    size_t num_nodes);

  void ExportCellResidualMoments();

  CellImportanceInfo
    GetCellImportanceInfo(const chi_mesh::Cell& cell, size_t g) const;

}; //MCPARTRA_RMC_SOURCE_A_H



#endif