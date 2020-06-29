#ifndef _mc_rmc3_source_h
#define _mc_rmc3_source_h

#include "../mc_base_source.h"

#include <ChiPhysics/chi_physics.h>

#include <ChiMath/chi_math.h>

//###################################################################
/**Residual source class.*/
class chi_montecarlon::ResidualSource3 : public chi_montecarlon::Source
{
public:
  chi_physics::FieldFunction* resid_ff;
private:
  struct CellSideGeometryData
  {
    double volume=0.0;
    double area=0.0;
    int    associated_face=-1;
    chi_mesh::Vector3 ref_point;
    std::vector<chi_mesh::Vector3> legs;
  };
  struct CellGeometryData
  {
    double total_volume=0.0;
    double total_area=0.0;
    std::vector<CellSideGeometryData> sides;
  };
  std::vector<CellGeometryData> cell_geometry_info;
public:
  std::vector<double> cell_avg_interior_rstar;
  std::vector<double> cell_avg_surface_rstar;
  std::vector<double> cell_volumes;
  std::vector<double> cell_IntVOmega_rstar;
  std::vector<double> cell_IntSOmega_rstar;

  double IntVOmega_rstar = 0.0;
  double IntSOmega_rstar = 0.0;
  double domain_volume = 0.0;
  double source_volume = 0.0;

  std::vector<double> interior_cdf;
  std::vector<double> surface_cdf;


  const bool sample_uniformly;
public:
  ResidualSource3(chi_physics::FieldFunction* in_resid_ff,
                  bool use_uniform_sampling=false);

  void Initialize(chi_mesh::MeshContinuum* ref_grid,
                  SpatialDiscretization_FV*   ref_fv_sdm,
                  chi_montecarlon::Solver* ref_solver);

  void BuildCellVolInfo(chi_mesh::MeshContinuum*  ref_grid,
                        SpatialDiscretization_FV* ref_fv_sdm);
  chi_mesh::Vector3 GetRandomPositionInCell(
    chi_montecarlon::RandomNumberGenerator& rng,
    const CellGeometryData& cell_info);
  chi_mesh::Vector3 GetRandomPositionOnCellSurface(
    chi_montecarlon::RandomNumberGenerator& rng,
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
    chi_montecarlon::RandomNumberGenerator& rng);
  chi_mesh::Vector3 RandomCosineLawDirection(
    chi_montecarlon::RandomNumberGenerator& rng,
    const chi_mesh::Vector3& normal);


  chi_montecarlon::Particle
  CreateParticle(chi_montecarlon::RandomNumberGenerator* rng);



};



#endif