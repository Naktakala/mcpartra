#include "mc_rmcA_source.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <ChiMath/SpatialDiscretization/FiniteVolume/fv.h>

#include <ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include <ChiPhysics/chi_physics.h>
#include "chi_log.h"
#include "SourceDrivenSolver/sdsolver.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

#include <tuple>

//###################################################################
/**Build cell volume information.*/
void mcpartra::ResidualSourceA::BuildCellVolInfo(
  chi_mesh::MeshContinuumPtr ref_grid,
  std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm)
{
//  for (auto cell_g_index : ref_grid->local_cell_glob_indices)
//  {
//    auto& cell = ref_grid->cells[cell_g_index];
  for (auto& cell : ref_grid->local_cells)
  {
    auto fv_view = ref_fv_sdm->MapFeView(cell.local_id);

    if (cell.Type() == chi_mesh::CellType::SLAB)
    {
      CellGeometryData cell_info;

      {
        const auto& v0 = ref_grid->vertices[cell.vertex_ids[0]];
        const auto& v1 = ref_grid->vertices[cell.vertex_ids[1]];

        auto v01 = v1 - v0;

        double V = 0.5*v01.Norm();

        CellSideGeometryData side_data;

        side_data.volume = V;
        side_data.area   = 0.5;
        side_data.associated_face = 0;
        side_data.ref_point = v0;
        side_data.legs.push_back(v01);

        cell_info.total_volume += V; //total volume
        cell_info.total_area += side_data.area;
        cell_info.sides.push_back(side_data);
      }
      {
        const auto& v0 = ref_grid->vertices[cell.vertex_ids[0]];
        const auto& v1 = ref_grid->vertices[cell.vertex_ids[1]];

        auto v10 = v0 - v1;

        double V = 0.5*v10.Norm();

        CellSideGeometryData side_data;

        side_data.volume = V;
        side_data.area   = 0.5;
        side_data.associated_face = 1;
        side_data.ref_point = v1;
        side_data.legs.push_back(v10);

        cell_info.total_volume += V; //total volume
        cell_info.total_area += side_data.area;
        cell_info.sides.push_back(side_data);
      }

      cell_geometry_info.push_back(cell_info);
    }
    else if (cell.Type() == chi_mesh::CellType::POLYGON)
    {
      CellGeometryData cell_info;
      int f=-1;
      for (auto& face : cell.faces)
      {
        ++f;
        const auto& v0 = ref_grid->vertices[face.vertex_ids[0]];
        const auto& v1 = ref_grid->vertices[face.vertex_ids[1]];
        auto& v2 = cell.centroid;

        auto v01 = v1-v0;
        auto v02 = v2-v0;

        double V = ((v01.x)*(v02.y) - (v02.x)*(v01.y))/2.0;

        CellSideGeometryData side_data;
        side_data.volume = V;
        side_data.area   = v01.Norm();
        side_data.associated_face = f;
        side_data.ref_point = v0;
        side_data.legs.push_back(v01);
        side_data.legs.push_back(v02);

        cell_info.total_volume += V;
        cell_info.total_area += side_data.area;
        cell_info.sides.push_back(side_data);
      }//for face

      cell_geometry_info.push_back(cell_info);
    }
    else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto& polyh_cell = (chi_mesh::CellPolyhedron&)cell;
      CellGeometryData cell_info;

      int f=-1;
      for (auto& face : cell.faces)
      {
        ++f;
        auto face_edges = polyh_cell.GetFaceEdges(f);

        for (auto& edge : face_edges)
        {
          const auto& v0 = ref_grid->vertices[edge[0]];
          const auto& v1 = ref_grid->vertices[edge[1]];
          auto& v2 = face.centroid;
          auto& v3 = cell.centroid;

          auto v01 = v1-v0;
          auto v02 = v2-v0;
          auto v03 = v3-v0;

          chi_mesh::Matrix3x3 J;
          J.SetColJVec(0,v01);
          J.SetColJVec(1,v02);
          J.SetColJVec(2,v03);

          double detJ = J.Det();

          double V = detJ/6.0;

          //Compute face portion area
          double a_00 = v01.Dot(v01);
          double a_01 = v01.Dot(v02);
          double a_10 = v02.Dot(v01);
          double a_11 = v02.Dot(v02);

          double A_f = 0.5*std::sqrt(a_00*a_11 - a_01*a_10);

          CellSideGeometryData side_data;
          side_data.volume = V;
          side_data.area   = A_f;
          side_data.associated_face = f;
          side_data.ref_point = v0;
          side_data.legs.push_back(v01);
          side_data.legs.push_back(v02);
          side_data.legs.push_back(v03);

          cell_info.total_volume += V;
          cell_info.total_area += A_f;
          cell_info.sides.push_back(side_data);
        }//for edge
      }//for face

      cell_geometry_info.push_back(cell_info);
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
/**Populates a material information data structure from a mat-id.*/
void mcpartra::ResidualSourceA::
  PopulateMaterialData(int mat_id, int group_g, MaterialData &mat_data)
{
  int  xs_prop_id     = ref_solver.matid_xs_map[mat_id];
  int  src_prop_id    = ref_solver.matid_q_map[mat_id];
  auto material = chi_physics_handler.material_stack[mat_id];
  auto xs = std::static_pointer_cast<chi_physics::TransportCrossSections>(
    material->properties[xs_prop_id]);

  double siga = xs->sigma_a[group_g];
  double Q    = 0.0;
  if (src_prop_id >= 0)
  {
    auto prop = material->properties[src_prop_id];
    auto q_prop =
      std::static_pointer_cast<chi_physics::IsotropicMultiGrpSource>(prop);
    Q = q_prop->source_value_g[group_g];
  }

  mat_data.siga = siga;
  mat_data.Q = Q;
}

//###################################################################
/**Samples the cell interior*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  GetRandomPositionInCell(
    chi_math::RandomNumberGenerator& rng,
    const mcpartra::ResidualSourceA::CellGeometryData &cell_info)
{
  chi_mesh::Vector3 position;
  double rn = rng.Rand();
  double cell_volume = cell_info.total_volume;
  double incremental_volume = 0.0;

  CellSideGeometryData const* ref_side_data = &cell_info.sides.back();
  for (auto& side_data : cell_info.sides)
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
    double u = rng.Rand();
    position = ref_side_data->ref_point +
               ref_side_data->legs[0]*u;
  }
  else if (num_legs == 2)
  {
    double u=rng.Rand();
    double v=rng.Rand();

    while ((u+v)>1.0)
    {u=rng.Rand(); v=rng.Rand();}

    position =
      ref_side_data->ref_point +
      ref_side_data->legs[0]*u +
      ref_side_data->legs[1]*v;
  }
  else if (num_legs == 3)
  {
    double u=rng.Rand();
    double v=rng.Rand();
    double w=rng.Rand();

    while ((u+v+w)>1.0)
    {u=rng.Rand(); v=rng.Rand(); w=rng.Rand();}

    position = ref_side_data->ref_point +
               ref_side_data->legs[0]*u +
               ref_side_data->legs[1]*v +
               ref_side_data->legs[2]*w;
  }


  return position;
}

//###################################################################
/**Samples the cell interior*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  GetRandomPositionOnCellSurface(
    chi_math::RandomNumberGenerator& rng,
    const mcpartra::ResidualSourceA::CellGeometryData &cell_info,
    const int face_mask,
    int* face_sampled)
{
  chi_mesh::Vector3 position;
  double rn = rng.Rand();
  double cell_area = cell_info.total_area;
  double incremental_area = 0.0;

  //======================================== Determine side
  CellSideGeometryData const* ref_side_data = &cell_info.sides.back();
  //=================================== Random face
  if (face_mask<0)
  {
    for (auto& side_data : cell_info.sides)
    {
      incremental_area += side_data.area;
      if (rn <= (incremental_area/cell_area))
      {
        ref_side_data = &side_data;
        break;
      }//if rn<cdf
    }//for side
    if (face_sampled != nullptr)
      *face_sampled = ref_side_data->associated_face;
  }
  //=================================== Specific face
  else
  {
    //============================ Build list of sides for face
    std::vector<CellSideGeometryData const*> face_sides;
    double face_area=0.0;
    for (const auto& side_data : cell_info.sides)
      if (side_data.associated_face == face_mask)
      {
        face_area += side_data.area;
        face_sides.push_back(&side_data);
      }

    //============================ Sample from list
    for (auto side_data : face_sides)
    {
      incremental_area += side_data->area;
      if (rn <= (incremental_area/face_area))
      {
        ref_side_data = side_data;
        break;
      }//if rn<cdf
    }//for side

    if (face_sampled != nullptr)
      *face_sampled = face_mask;
  }

  //======================================== Sample face
  auto num_legs = ref_side_data->legs.size();

  if (num_legs == 1)
  {
    position = ref_side_data->ref_point;
  }
  else if (num_legs == 2)
  {
    double u=rng.Rand();

    position =
      ref_side_data->ref_point +
      ref_side_data->legs[0]*u;
  }
  else if (num_legs == 3)
  {
    double u=rng.Rand();
    double v=rng.Rand();

    while ((u+v)>1.0)
    {u=rng.Rand(); v=rng.Rand();}

    position =
      ref_side_data->ref_point +
      ref_side_data->legs[0]*u +
      ref_side_data->legs[1]*v;
  }


  return position;
}

//###################################################################
/**Obtains a field function interpolant of the flux.*/
double mcpartra::ResidualSourceA::
GetResidualFFPhi(std::vector<double> &N_in,
                 int dofs,
                 int cell_local_id,
                 int egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  double phi = 0.0;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = pwl_sdm->MapDOFLocal(cell,dof,uk_man,0,egrp);

    phi += (*resid_ff->field_vector_local)[ir]*N_in[dof];
  }//for dof

  return phi;
}

//###################################################################
/**Obtains a field function interpolant of the flux-gradient.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
GetResidualFFGradPhi(std::vector<chi_mesh::Vector3>& Grad_in,
                     int dofs,
                     int cell_local_id,
                     int egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  chi_mesh::Vector3 gradphi;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = pwl_sdm->MapDOFLocal(cell,dof,uk_man,0,egrp);

    gradphi = gradphi + (Grad_in[dof]*(*resid_ff->field_vector_local)[ir]);
  }//for dof

  return gradphi;
}

//###################################################################
/**Gets a true random direction.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  RandomDirection(
    chi_math::RandomNumberGenerator& rng)
{
  double costheta = 2.0*rng.Rand()-1.0;
  double theta    = acos(costheta);
  double varphi   = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 omega;
  omega.x = sin(theta)*cos(varphi);
  omega.y = sin(theta)*sin(varphi);
  omega.z = cos(theta);

  return omega;
}

//###################################################################
/**Gets a cosine law random direction relative to a normal.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  RandomCosineLawDirection(
    chi_math::RandomNumberGenerator& rng,
    const chi_mesh::Vector3& normal)
{
  //Build rotation matrix
  chi_mesh::Matrix3x3 R;

  chi_mesh::Vector3 khat(0.0,0.0,1.0);

  if      (normal.Dot(khat) >  0.9999999)
    R.SetDiagonalVec(1.0,1.0,1.0);
  else if (normal.Dot(khat) < -0.9999999)
    R.SetDiagonalVec(1.0,1.0,-1.0);
  else
  {
    chi_mesh::Vector3 binorm = khat.Cross(normal);
    binorm = binorm/binorm.Norm();

    chi_mesh::Vector3 tangent = binorm.Cross(normal);
    tangent = tangent/tangent.Norm();

    R.SetColJVec(0,tangent);
    R.SetColJVec(1,binorm);
    R.SetColJVec(2,normal);
  }

  //Sample direction
  double costheta = rng.Rand();     //Sample half-range only
  double theta    = acos(sqrt(costheta));
  double varphi   = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  return R*ref_dir;
}