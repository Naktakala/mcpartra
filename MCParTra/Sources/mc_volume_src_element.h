#ifndef MCPARTA_VOLUME_SOURCE_ELEMENT_H
#define MCPARTA_VOLUME_SOURCE_ELEMENT_H

#include "../mcpartra.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include <cstdint>
#include <utility>

namespace mcpartra
{
//###################################################################
/**A basic element that can be used to generate particles.*/
class VolumeSourceElement
{
private:
  const uint64_t parent_cell_local_id;
  const uint64_t parent_cell_global_id;

  const chi_mesh::Vector3 ref_point;
  const std::vector<chi_mesh::Vector3> geom_legs;

  double volume = 0.0;

public:
  VolumeSourceElement(uint64_t in_cell_local_id,
                      uint64_t in_cell_global_id,
                      const chi_mesh::Vector3& in_ref_point,
                      std::vector<chi_mesh::Vector3> in_geom_legs);

  chi_mesh::Vector3 SampleRandomPosition(chi_math::RandomNumberGenerator& rng);

  double Volume() const {return volume;}
  uint64_t ParentCellLocalID() const {return parent_cell_local_id;}
  uint64_t ParentCellGlobalID() const {return parent_cell_global_id;}
  size_t TypeIndex() const {return geom_legs.size();}
};

//###################################################################
/**Make a list of VolumeSourceElement-type elements for a given cell.*/
std::vector<VolumeSourceElement>
  GetCellVolumeSourceElements(chi_mesh::Cell& cell,
                              chi_mesh::MeshContinuumPtr& grid);

}//namespace mcpartra

#endif //MCPARTA_VOLUME_SOURCE_ELEMENT_H