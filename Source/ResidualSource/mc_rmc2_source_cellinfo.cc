#include "mc_rmc2_source.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include "../../Solver/solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Build cell volume information.*/
void chi_montecarlon::ResidualSource2::BuildCellVolInfo(
  chi_mesh::MeshContinuum* ref_grid,
  SpatialDiscretization_FV* ref_fv_sdm)
{
  for (auto cell_g_index : ref_grid->local_cell_glob_indices)
  {
    auto cell = ref_grid->cells[cell_g_index];
    auto fv_view = ref_fv_sdm->MapFeView(cell->local_id);

    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      CellSideInfo cell_side_info;

      auto v0 = *ref_grid->vertices[cell->vertex_ids[0]];
      auto v1 = *ref_grid->vertices[cell->vertex_ids[1]];

      auto v01 = v1 - v0;

      double V = fv_view->volume;

      CellSideData side_data;

      side_data.volume = V;
      side_data.ref_point = v0;
      side_data.legs.push_back(v01);

      cell_side_info.first = V; //total volume
      cell_side_info.second.push_back(side_data);

      cell_vol_info.push_back(cell_side_info);
    }
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      CellSideInfo cell_side_info(0.0,std::vector<CellSideData>());
      for (auto& face : cell->faces)
      {
        auto  v0 = *ref_grid->vertices[face.vertex_ids[0]];
        auto  v1 = *ref_grid->vertices[face.vertex_ids[1]];
        auto& v2 = cell->centroid;

        auto v01 = v1-v0;
        auto v02 = v2-v0;

        double V = ((v01.x)*(v02.y) - (v02.x)*(v01.y))/2.0;

        CellSideData side_data;
        side_data.volume = V;
        side_data.ref_point = v0;
        side_data.legs.push_back(v01);
        side_data.legs.push_back(v02);

        cell_side_info.first += V;
        cell_side_info.second.push_back(side_data);
      }//for face

      cell_vol_info.push_back(cell_side_info);
    }
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
      size_t f=0;
      CellSideInfo cell_side_info(0.0,std::vector<CellSideData>());
      for (auto& face : cell->faces)
      {
        auto face_edges = polyh_cell->GetFaceEdges(f);

        for (auto& edge : face_edges)
        {
          auto  v0 = *ref_grid->vertices[edge[0]];
          auto  v1 = *ref_grid->vertices[edge[1]];
          auto& v2 = face.centroid;
          auto& v3 = cell->centroid;

          auto v01 = v1-v0;
          auto v02 = v2-v0;
          auto v03 = v3-v0;

          chi_mesh::Matrix3x3 J;
          J.SetColJVec(0,v01);
          J.SetColJVec(1,v02);
          J.SetColJVec(2,v03);

          double detJ = J.Det();

          double V = detJ/6.0;

          CellSideData side_data;
          side_data.volume = V;
          side_data.ref_point = v0;
          side_data.legs.push_back(v01);
          side_data.legs.push_back(v02);
          side_data.legs.push_back(v03);

          cell_side_info.first += V;
          cell_side_info.second.push_back(side_data);
        }//for edge
        ++f;
      }//for face

      cell_vol_info.push_back(cell_side_info);
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encountered in "
        << "chi_montecarlon::MaterialSource:: Initialize.";
      exit(EXIT_FAILURE);
    }
  }//for local cell

}

//###################################################################
/**Samples the cell interior*/
chi_mesh::Vector3 chi_montecarlon::ResidualSource2::
  GetRandomPositionInCell(
    chi_math::RandomNumberGenerator *rng,
    chi_montecarlon::ResidualSource2::CellSideInfo &cell_side_info)
{
  chi_mesh::Vector3 position;
  double rn = rng->Rand();
  double cell_volume = cell_side_info.first;
  double incremental_volume = 0.0;

  CellSideData* ref_side_data = &cell_side_info.second.back();
  for (auto& side_data : cell_side_info.second)
  {
    incremental_volume += side_data.volume;
    if (rn <= (incremental_volume/cell_volume))
    {
      ref_side_data = &side_data;
      break;
    }//if rn<cdf
  }//for side

  auto num_legs = ref_side_data->legs.size();

  if (num_legs == 1)
  {
    double u = rng->Rand();
    position = ref_side_data->ref_point +
               ref_side_data->legs[0]*u;
  }
  else if (num_legs == 2)
  {
    double u=rng->Rand();
    double v=rng->Rand();

    while ((u+v)>1.0)
    {u=rng->Rand(); v=rng->Rand();}

    position =
      ref_side_data->ref_point +
      ref_side_data->legs[0]*u +
      ref_side_data->legs[1]*v;
  }
  else if (num_legs == 3)
  {
    double u=rng->Rand();
    double v=rng->Rand();
    double w=rng->Rand();

    while ((u+v+w)>1.0)
    {u=rng->Rand(); v=rng->Rand(); w=rng->Rand();}

    position = ref_side_data->ref_point +
               ref_side_data->legs[0]*u +
               ref_side_data->legs[1]*v +
               ref_side_data->legs[2]*w;
  }


  return position;
}