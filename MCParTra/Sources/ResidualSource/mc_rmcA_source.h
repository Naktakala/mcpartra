#ifndef MCPARTRA_RMC_SOURCE_A_H
#define MCPARTRA_RMC_SOURCE_A_H

#include "../mc_base_source.h"

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
  /**Structure to hold a cell's constitute primitive info.
   * - For a slab this is just a duality to get points on each face.
   * - For a polygon this is the constituent triangles.
   * - For a polyhedron this is the constituent tetrahedrons.*/
  struct CellSideGeometryData
  {
    double volume=0.0;
    double area=0.0;
    int    associated_face=-1;
    chi_mesh::Vector3 ref_point;
    std::vector<chi_mesh::Vector3> legs;
  };
  /**Structure to hold all of the constituents.*/
  struct CellGeometryData
  {
    double total_volume=0.0;
    double total_area=0.0;
    std::vector<CellSideGeometryData> sides;
  };
  std::vector<CellGeometryData> cell_geometry_info;

  /**Structure to store cell-face pairs.*/
  struct RCellFace
  {
    int cell_local_id=-1;
    int ass_face=-1;
    double average_rstar=0.0;
    double maximum=-1.0e32;
    double minimum= 1.0e32;
    double area=0.0;
  };
  std::vector<RCellFace> r_abs_cellk_facef_surface_average;

public:
  std::vector<double> r_abs_cellk_interior_average;
  std::vector<double> r_cellk_interior_max;
  std::vector<double> r_cellk_interior_min;
  std::vector<double> R_abs_cellk_interior;

  double R_abs_localdomain_interior = 0.0;
  double R_abs_localdomain_surface = 0.0;

  double R_abs_globaldomain_interior = 0.0;
  double R_abs_globaldomain_surface = 0.0;

  double relative_weight = 1.0;

  double domain_volume = 0.0;
  double source_volume = 0.0;

  std::vector<double> domain_cdf;
  std::vector<double> surface_cdf;


  const bool sample_uniformly;
public:
  explicit
  ResidualSourceA(mcpartra::SourceDrivenSolver& solver,
                  std::shared_ptr<chi_physics::FieldFunction>& in_resid_ff,
                  bool use_uniform_sampling=false) :
    SourceBase(SourceType::RESIDUAL_TYPE_A, solver),
    resid_ff(in_resid_ff),
    sample_uniformly(use_uniform_sampling)
  {}

  //a
  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm) override;

  //b
  void BuildCellVolInfo(chi_mesh::MeshContinuumPtr  ref_grid,
                        std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm);

  void PopulateMaterialData(int mat_id, int group_g,
                            MaterialData& mat_data);

  static chi_mesh::Vector3 GetRandomPositionInCell(
    chi_math::RandomNumberGenerator& rng,
    const CellGeometryData& cell_info);

  static chi_mesh::Vector3 GetRandomPositionOnCellSurface(
    chi_math::RandomNumberGenerator& rng,
    const CellGeometryData& cell_info,
    const int face_mask=-1,
    int* face_sampled= nullptr);

  double GetResidualFFPhi(std::vector<double> &N_in,
                          int dofs,
                          int cell_local_id,
                          int egrp);

  chi_mesh::Vector3 GetResidualFFGradPhi(
                          std::vector<chi_mesh::Vector3>& Grad_in,
                          int dofs,
                          int cell_local_id,
                          int egrp);

  chi_mesh::Vector3 RandomDirection(
    chi_math::RandomNumberGenerator& rng);

  chi_mesh::Vector3 RandomCosineLawDirection(
    chi_math::RandomNumberGenerator& rng,
    const chi_mesh::Vector3& normal);

  //c
  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

  //d
  void RemoveFFDiscontinuities();

}; //MCPARTRA_RMC_SOURCE_A_H



#endif