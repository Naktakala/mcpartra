#include "mc_volume_src_element.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Constructor.*/
mcpartra::VolumeSourceElement::
  VolumeSourceElement(uint64_t in_cell_local_id,
                      uint64_t in_cell_global_id,
                      const chi_mesh::Vector3& in_ref_point,
                      std::vector<chi_mesh::Vector3> in_geom_legs) :
    parent_cell_local_id(in_cell_local_id),
    parent_cell_global_id(in_cell_global_id),
    ref_point(in_ref_point),
    geom_legs(std::move(in_geom_legs))
{
  if (geom_legs.size() == 1) //Line
  {
    volume = geom_legs.back().Norm();
  }
  else if (geom_legs.size() == 2) //Triangle
  {
    const auto& v01 = geom_legs[0];
    const auto& v02 = geom_legs[1];

    volume = ((v01.x)*(v02.y) - (v02.x)*(v01.y))/2.0;
  }
  else if (geom_legs.size() == 3) //Tetrahedron
  {
    const auto& v01 = geom_legs[0];
    const auto& v02 = geom_legs[1];
    const auto& v03 = geom_legs[2];

    chi_mesh::Matrix3x3 J;
    J.SetColJVec(0,v01);
    J.SetColJVec(1,v02);
    J.SetColJVec(2,v03);

    double detJ = J.Det();

    volume = detJ/6.0;
  }
}

//###################################################################
/**Sample a random position within an elment.*/
chi_mesh::Vector3 mcpartra::VolumeSourceElement::
  SampleRandomPosition(chi_math::RandomNumberGenerator& rng)
{
  if (geom_legs.size() == 1) //Line
  {
    double u = rng.Rand();
    return ref_point + geom_legs[0]*u;
  }
  else if (geom_legs.size() == 2) //Polygon
  {
    double u = rng.Rand();
    double v = rng.Rand();

    while ((u+v)>1.0)
    {u=rng.Rand(); v=rng.Rand();}

    return ref_point + geom_legs[0]*u + geom_legs[1]*v;
  }
  else if (geom_legs.size() == 3) //Polyhedron
  {
    double u=rng.Rand();
    double v=rng.Rand();
    double w=rng.Rand();

    while ((u+v+w)>1.0)
    {u=rng.Rand(); v=rng.Rand(); w=rng.Rand();}

    return ref_point + geom_legs[0]*u + geom_legs[1]*v + geom_legs[2]*w;
  }

  return ref_point;
}

//###################################################################
/**Makes VolumeSourceElements for a given cell.*/
std::vector<mcpartra::VolumeSourceElement> mcpartra::
  GetCellVolumeSourceElements(chi_mesh::Cell &cell,
                              chi_mesh::MeshContinuumPtr& grid)
{
  std::vector<VolumeSourceElement> elements;

  if (cell.Type() == chi_mesh::CellType::SLAB)
  {
    const auto& v0 = grid->vertices[cell.vertex_ids[0]];
    const auto& v1 = grid->vertices[cell.vertex_ids[1]];

    auto v01 = v1 - v0;

    std::vector<chi_mesh::Vector3> legs{v01};

    elements.emplace_back(cell.local_id, cell.global_id, v0, legs);
  }
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    for (auto& face : cell.faces)
    {
      const auto& v0 = grid->vertices[face.vertex_ids[0]];
      const auto& v1 = grid->vertices[face.vertex_ids[1]];
      auto& v2 = cell.centroid;

      auto v01 = v1-v0;
      auto v02 = v2-v0;

      std::vector<chi_mesh::Vector3> legs{v01,v02};

      elements.emplace_back(cell.local_id, cell.global_id, v0, legs);
    }//for face
  }
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    for (auto& face : cell.faces)
    {
      size_t num_face_verts = face.vertex_ids.size();
      for (int fv=0; fv<num_face_verts; ++fv)
      {
        int fvp1 = (fv<(num_face_verts-1))? fv+1 : 0;
        const auto& v0 = grid->vertices[fv];
        const auto& v1 = grid->vertices[fvp1];
        auto& v2 = face.centroid;
        auto& v3 = cell.centroid;

        auto v01 = v1-v0;
        auto v02 = v2-v0;
        auto v03 = v3-v0;

        std::vector<chi_mesh::Vector3> legs{v01,v02,v03};

        elements.emplace_back(cell.local_id, cell.global_id, v0, legs);
      }//for edge
    }//for face
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Unsupported cell type encountered in "
      << "mcpartra::GetCellVolumeSourceElements.";
    exit(EXIT_FAILURE);
  }

  return elements;
}