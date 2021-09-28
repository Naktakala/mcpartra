#ifndef MCPARTRA_SURFACE_SRC_ELEMENT_H
#define MCPARTRA_SURFACE_SRC_ELEMENT_H

#include "../mcpartra.h"

#include "ChiMesh/chi_mesh.h"
#include "ChiMath/RandomNumberGeneration/random_number_generator.h"

#include "ChiMesh/Cell/cell.h"

#include <cstdint>
#include <utility>

namespace mcpartra
{
//###################################################################
/**A basic element that can be used to generate particles.*/
class SurfaceSourceElement
{
private:
  const uint64_t parent_cell_local_id;
  const uint64_t parent_cell_global_id;

  const chi_mesh::Vector3 ref_point;
  const std::vector<chi_mesh::Vector3> geom_legs;
  chi_mesh::Vector3 normal;

  double area = 0.0;

public:
  SurfaceSourceElement(uint64_t in_cell_local_id,
                       uint64_t in_cell_global_id,
                       const chi_mesh::Vector3& in_ref_point,
                       std::vector<chi_mesh::Vector3> in_geom_legs,
                       const chi_mesh::Vector3& in_normal);

  chi_mesh::Vector3 SampleRandomPosition(chi_math::RandomNumberGenerator& rng) const;

  double Area() const {return area;}
  chi_mesh::Vector3 Normal() const {return normal;}
  uint64_t ParentCellLocalID() const {return parent_cell_local_id;}
  uint64_t ParentCellGlobalID() const {return parent_cell_global_id;}
  size_t TypeIndex() const {return geom_legs.size();}
};

//###################################################################
/**Make a list of VolumeSourceElement-type elements for a given cell
 * comprising the surface elements.*/
std::vector<SurfaceSourceElement>
GetCellSurfaceSourceElements(const chi_mesh::Cell& cell,
                             const chi_mesh::CellFace& face,
                             const chi_mesh::MeshContinuumPtr& grid);
}//namespace mcpartra
#endif //MCPARTRA_SURFACE_SRC_ELEMENT_H