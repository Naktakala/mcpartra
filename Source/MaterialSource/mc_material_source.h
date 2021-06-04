#ifndef _montecarlon_material_source_h
#define _montecarlon_material_source_h

#include "../mc_base_source.h"

//###################################################################
/**Material source class.*/
class chi_montecarlon::MaterialSource : public chi_montecarlon::Source
{
private:
  struct CellSourceElement
  {
    const int cell_local_index;
    const int cell_global_index;
    const int group;
    const double product_V_Q;
    const double source_value;

    const chi_mesh::Vector3 ref_point;
    const std::vector<chi_mesh::Vector3> geom_legs;

    CellSourceElement(int in_cell_local_index,
                      int in_cell_global_index,
                      int in_group,
                      double in_product_V_Q,
                      double in_source_value,
                      const chi_mesh::Vector3& in_source_geom_ref_point,
                      std::vector<chi_mesh::Vector3>& in_source_geom_legs) :
      cell_local_index(in_cell_local_index),
      cell_global_index(in_cell_global_index),
      group(in_group),
      product_V_Q(in_product_V_Q),
      source_value(in_source_value),
      ref_point(in_source_geom_ref_point),
      geom_legs(in_source_geom_legs)
    { }
  };
  std::vector<double> IntV_Q_g;
  std::vector<double> group_cdf;
  std::vector<std::vector<CellSourceElement>> group_sources;
  std::vector<std::vector<double>> group_element_cdf;



public:
  MaterialSource() {};

  void Initialize(chi_mesh::MeshContinuumPtr ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm,
                  chi_montecarlon::Solver* ref_solver) override;

  chi_montecarlon::Particle
  CreateParticle(chi_math::RandomNumberGenerator* rng) override;

  double GetParallelRelativeSourceWeight() override;
};

#endif