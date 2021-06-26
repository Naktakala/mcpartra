#ifndef MCPARTRA_MATERIAL_SOURCE_H
#define MCPARTRA_MATERIAL_SOURCE_H

#include "../mc_base_source.h"
#include "../mc_volume_src_element.h"

//###################################################################
/**Material source class.*/
class mcpartra::MaterialSource : public mcpartra::SourceBase
{
private:
  struct CellSourceElement
  {
    const int cell_local_id;
    const int cell_global_id;
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
      cell_local_id(in_cell_local_index),
      cell_global_id(in_cell_global_index),
      group(in_group),
      product_V_Q(in_product_V_Q),
      source_value(in_source_value),
      ref_point(in_source_geom_ref_point),
      geom_legs(in_source_geom_legs)
    { }
  };
  std::vector<double> IntV_Q_g;
  std::vector<double> group_cdf;
  std::map<uint64_t, std::vector<VolumeSourceElement>> cell_elements;

  typedef std::pair<double, VolumeSourceElement*> ElementSrc;
  typedef std::vector<ElementSrc> GrpSrc;

  std::vector<GrpSrc> group_sources;
  std::vector<std::vector<double>> group_element_cdf;



public:
  explicit
  MaterialSource(mcpartra::Solver& solver) :
    SourceBase(SourceType::MATERIAL_SRC, solver)
  {};

  void Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
                  std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm) override;

  mcpartra::Particle
  CreateParticle(chi_math::RandomNumberGenerator& rng) override;

//  double GetParallelRelativeSourceWeight() override;
};

#endif //MCPARTRA_MATERIAL_SOURCE_H