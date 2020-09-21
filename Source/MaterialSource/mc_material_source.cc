#include "mc_material_source.h"

#include "../../Solver/solver_montecarlon.h"

#include "ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Initialize material source.*/
void chi_montecarlon::MaterialSource::
  Initialize(chi_mesh::MeshContinuum *ref_grid,
             SpatialDiscretization_FV *ref_fv_sdm,
             chi_montecarlon::Solver *ref_solver)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  this->ref_solver = ref_solver;

  //============================================= List cells with mat-sources
  std::vector<chi_mesh::Cell*> mat_src_cells;
  for (auto& cell : ref_grid->local_cells)
  {
    auto mat = chi_physics_handler.material_stack[cell.material_id];
    for (auto prop : mat->properties)
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        mat_src_cells.push_back(&cell);
        break;
      }
  }

  //============================================= Determine group-wise
  //                                              source weights
  IntV_Q_g.clear();
  IntV_Q_g.resize(ref_solver->num_grps,0.0);
  group_sources.resize(ref_solver->num_grps);
  for (auto cell : mat_src_cells)
  {
    auto fv_view = ref_fv_sdm->MapFeView(cell->local_id);

    auto mat = chi_physics_handler.material_stack[cell->material_id];

    for (auto prop : mat->properties)
    {
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = (chi_physics::IsotropicMultiGrpSource*)prop;
        for (size_t g=0; g<src->source_value_g.size(); ++g)
        {
          IntV_Q_g[g] += fv_view->volume*src->source_value_g[g];

          if (cell->Type() == chi_mesh::CellType::SLAB)
          {
            auto v0 = *ref_grid->vertices[cell->vertex_ids[0]];
            auto v1 = *ref_grid->vertices[cell->vertex_ids[1]];

            auto v01 = v1 - v0;

            double V = fv_view->volume;

            std::vector<chi_mesh::Vector3> legs;
            legs.push_back(v01);

            if (src->source_value_g[g]>1.0e-12)
              group_sources[g].emplace_back(cell->local_id,
                                            cell->global_id,
                                            g,
                                            V*src->source_value_g[g],
                                            src->source_value_g[g],
                                            v0,
                                            legs);
          }
          else if (cell->Type() == chi_mesh::CellType::POLYGON)
          {
            for (auto& face : cell->faces)
            {
              auto  v0 = *ref_grid->vertices[face.vertex_ids[0]];
              auto  v1 = *ref_grid->vertices[face.vertex_ids[1]];
              auto& v2 = cell->centroid;

              auto v01 = v1-v0;
              auto v02 = v2-v0;

              double V = ((v01.x)*(v02.y) - (v02.x)*(v01.y))/2.0;

              std::vector<chi_mesh::Vector3> legs;
              legs.push_back(v01);
              legs.push_back(v02);

              if (src->source_value_g[g]>1.0e-12)
                group_sources[g].emplace_back(cell->local_id,
                                              cell->global_id,
                                              g,
                                              V*src->source_value_g[g],
                                              src->source_value_g[g],
                                              v0,
                                              legs);
            }//for face
          }
          else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
          {
            auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
            size_t f=0;
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

                std::vector<chi_mesh::Vector3> legs;
                legs.push_back(v01);
                legs.push_back(v02);
                legs.push_back(v03);

                if (src->source_value_g[g]>1.0e-12)
                  group_sources[g].emplace_back(cell->local_id,
                                                cell->global_id,
                                                g,
                                                V*src->source_value_g[g],
                                                src->source_value_g[g],
                                                v0,
                                                legs);
              }//for edge
              ++f;
            }//for face
          }
          else
          {
            chi_log.Log(LOG_ALLERROR)
              << "Unsupported cell type encountered in "
              << "chi_montecarlon::MaterialSource:: Initialize.";
            exit(EXIT_FAILURE);
          }
          break;
        }
      }//if src prop
    }
  }//for cell

  // The next two steps build two levels of CDFs. The first is group-wise,
  // then the second is the elements within a group.

  //============================================= Construct groupwise cdf
  double IntV_Q_total = 0.0;
  for (auto val : IntV_Q_g)
    IntV_Q_total += val;

  group_cdf.clear();
  group_cdf.resize(ref_solver->num_grps,0.0);
  double running_total = 0.0;
  for (size_t g=0; g<ref_solver->num_grps; ++g)
  {
    running_total += IntV_Q_g[g];
    group_cdf[g] = running_total/IntV_Q_total;
    chi_log.Log(LOG_0) << "Group cdf " << g << " " << group_cdf[g];
  }

  //============================================= Within each group, construct
  //                                              source element cdf
  group_element_cdf.clear();
  group_element_cdf.resize(ref_solver->num_grps);
  int g=0;
  for (auto& group_source : group_sources)
  {
    double group_total = 0.0;
    for (auto& src_element : group_source)
      group_total += src_element.product_V_Q;

    group_element_cdf[g].resize(group_source.size(),0.0);
    double elem_running_total=0.0;
    int elem=0;
    for (auto& src_element : group_source)
    {
      elem_running_total += src_element.product_V_Q;
      group_element_cdf[g][elem] = elem_running_total/group_total;

      ++elem;
    }
    ++g;
  }

}

//###################################################################
/**Creates a source particle.*/
chi_montecarlon::Particle chi_montecarlon::MaterialSource::
  CreateParticle(chi_math::RandomNumberGenerator *rng)
{
  Particle new_particle;

  //======================================== Sample group
  int g = std::lower_bound(
            group_cdf.begin(),
            group_cdf.end(),
            rng->Rand()) - group_cdf.begin();

  new_particle.egrp = g;

  if (group_element_cdf[g].empty())
  {
    new_particle.alive = false;
    return new_particle;
  }

  //======================================== Sample element
  int elem = std::lower_bound(
               group_element_cdf[g].begin(),
               group_element_cdf[g].end(),
               rng->Rand()) - group_element_cdf[g].begin();

  auto& src_element = group_sources[g][elem];

  //======================================== Sample position
  if      (src_element.geom_legs.size() == 1) //Slab
  {
    double u = rng->Rand();
    new_particle.pos = src_element.ref_point +
                       src_element.geom_legs[0]*u;
  }
  else if (src_element.geom_legs.size() == 2) //Triangle
  {
    double u=rng->Rand();
    double v=rng->Rand();

    while ((u+v)>1.0)
    {u=rng->Rand(); v=rng->Rand();}

    new_particle.pos = src_element.ref_point +
                       src_element.geom_legs[0]*u +
                       src_element.geom_legs[1]*v;
  }
  else if (src_element.geom_legs.size() == 3) //Tetrahedron
  {
    double u=rng->Rand();
    double v=rng->Rand();
    double w=rng->Rand();

    while ((u+v+w)>1.0)
    {u=rng->Rand(); v=rng->Rand(); w=rng->Rand();}

    new_particle.pos = src_element.ref_point +
                       src_element.geom_legs[0]*u +
                       src_element.geom_legs[1]*v +
                       src_element.geom_legs[2]*w;
  }

  //======================================== Sample direction
  double costheta = 2.0*rng->Rand() - 1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;

  //======================================== Determine weight
  new_particle.w = 1.0;

  new_particle.cur_cell_local_id  = src_element.cell_local_index;
  new_particle.cur_cell_global_id = src_element.cell_global_index;

  if (ref_solver->uncollided_only)
    new_particle.ray_trace_method = chi_montecarlon::Solver::RayTraceMethod::UNCOLLIDED;

  return new_particle;
}

//###################################################################
/**Gets the relative source strength accross all processors.*/
double chi_montecarlon::MaterialSource::GetParallelRelativeSourceWeight()
{
  double local_total_source_weight = 0.0;

  for (auto& grp_src : group_sources)
    for (auto& src_element : grp_src)
      local_total_source_weight += src_element.product_V_Q;

  double global_total_source_weight = 0.0;
  MPI_Allreduce(&local_total_source_weight,  //sendbuf
                &global_total_source_weight, //recvbuf
                1,                           //recvcount
                MPI_DOUBLE,                  //datatype
                MPI_SUM,                     //operation
                MPI_COMM_WORLD);             //communicator

  return local_total_source_weight/global_total_source_weight;
}