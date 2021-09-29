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
    double            Rstar_absolute=0.0;
  };
  typedef std::vector<RCellInterior> VecRCellInterior;
  std::vector<VecRCellInterior> residual_info_cell_interiors;

  /**Structure to store cell-face pairs.*/
  struct RCellFace
  {
    uint64_t cell_local_id=0;
    unsigned int ass_face=0;
    double maximum_rstar_absolute=-1.0e32;
    double Rstar_absolute=0.0;
  };
  typedef std::vector<RCellFace> VecRCellFace;
  std::vector<VecRCellFace> residual_info_cell_bndry_faces;

  /**Structure to store cell-face pairs.*/
  struct RCellFace2 : public RCellInterior
  {
    unsigned int ass_face=0;
    RCellFace2() {type = RessidualInfoType::Face;}
  };

private://PDFs
  std::vector<std::unique_ptr<RCellInterior>> residual_info_elements;

private:
  std::vector<double> R_abs_cellk_interior; ///< For field function

private: //CDFs
  std::vector<double> R_abs_domain_interior; ///< Per group
  std::vector<double> R_abs_domain_surface;  ///< Per group
  std::vector<double> R_abs_domain_total; ///< Per group

  std::vector<std::vector<double>> R_abs_domain_interior_cdf; ///< Per group then cell
  std::vector<std::vector<double>> R_abs_domain_surface_cdf; ///< Per group then face

  std::vector<std::vector<double>> interior_vs_surface_cdf; ///< Per group then exterior or surface

  std::vector<double> group_cdf; ///< Per group
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
  void PopulateMaterialData(int mat_id, size_t group_g,
                            MaterialData& mat_data);

  double GetResidualFFPhi(std::vector<double> &N_in,
                          size_t dofs,
                          uint64_t cell_local_id,
                          size_t egrp);

  chi_mesh::Vector3 GetResidualFFGradPhi(
                          std::vector<chi_mesh::Vector3>& Grad_in,
                          size_t dofs,
                          uint64_t cell_local_id,
                          size_t egrp);

  //c
  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

  //d
  void RemoveFFDiscontinuities();

}; //MCPARTRA_RMC_SOURCE_A_H



#endif