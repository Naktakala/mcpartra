#include "mc_material_source.h"

#include "../../Solver/solver_montecarlon.h"

#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Initialize material source.*/
void mcpartra::MaterialSource::
  Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm)
{
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;

  //============================================= List cells with mat-sources
  std::map<uint64_t, chi_mesh::Cell*> mat_src_cells;
  for (auto& cell : ref_grid->local_cells)
  {
    auto mat = chi_physics_handler.material_stack[cell.material_id];
    for (auto& prop : mat->properties)
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        mat_src_cells.insert(std::make_pair(cell.local_id,&cell));
        cell_elements.insert(
          std::make_pair(cell.local_id, GetCellVolumeSourceElements(cell,grid)));
        break;
      }
  }

  //============================================= Determine group-wise
  //                                              source weights
  typedef chi_physics::IsotropicMultiGrpSource IsoMGSrc;
  IntV_Q_g.clear();
  IntV_Q_g.resize(ref_solver.num_grps,0.0);
  group_sources.resize(ref_solver.num_grps);
  for (auto& cell_index_ptr_pair : mat_src_cells)
  {
    auto& cell = *cell_index_ptr_pair.second;
    auto fv_view = ref_fv_sdm->MapFeView(cell.local_id);

    auto mat = chi_physics_handler.material_stack[cell.material_id];

    for (auto& prop : mat->properties)
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto src = std::static_pointer_cast<IsoMGSrc>(prop);
        for (size_t g=0; g<src->source_value_g.size(); ++g)
        {
          double Q_g = src->source_value_g[g];
          IntV_Q_g[g] += fv_view->volume*Q_g;

          const auto& v_elements = cell_elements[cell.local_id];
          for (auto& element : v_elements)
            group_sources[g].emplace_back(Q_g * element.Volume(), element);
        }//for g
      }//if src prop
  }//for cell

  // The next two steps build two levels of CDFs. The first is group-wise,
  // then the second is the elements within a group.

  //============================================= Construct groupwise cdf
  double IntV_Q_total_local = 0.0;
  for (auto val : IntV_Q_g)
    IntV_Q_total_local += val;

  double IntV_Q_total_globl = 0.0;
  MPI_Allreduce(&IntV_Q_total_local,     // sendbuf
                &IntV_Q_total_globl,     // recvbuf
                1, MPI_DOUBLE,           // count + datatype
                MPI_SUM,                 // operation
                MPI_COMM_WORLD);         // communicator

  SetSourceRates(IntV_Q_total_local,IntV_Q_total_globl);

  group_cdf.clear();
  group_cdf.resize(ref_solver.num_grps,0.0);
  double running_total = 0.0;
  for (size_t g=0; g<ref_solver.num_grps; ++g)
  {
    running_total += IntV_Q_g[g];
    group_cdf[g] = running_total / IntV_Q_total_local;
  }

  //============================================= Within each group, construct
  //                                              source element cdf
  group_element_cdf.clear();
  group_element_cdf.resize(ref_solver.num_grps);
  int g=0;
  for (auto& group_source : group_sources)
  {
    group_element_cdf[g].resize(group_source.size(),0.0);
    double elem_running_total=0.0;
    size_t elem=0;
    for (auto& src_element : group_source)
    {
      elem_running_total += src_element.first;
      group_element_cdf[g][elem] = elem_running_total/IntV_Q_g[g];

      ++elem;
    }
    ++g;
  }

}


//###################################################################
/**Creates a source particle.*/
mcpartra::Particle mcpartra::MaterialSource::
  CreateParticle(chi_math::RandomNumberGenerator *rng)
{
  Particle new_particle;

  //======================================== Sample group
  size_t g = std::lower_bound(
            group_cdf.begin(),
            group_cdf.end(),
            rng->Rand()) - group_cdf.begin();

  if (g>group_cdf.size())
    return Particle::MakeDeadParticle();

  new_particle.egrp = static_cast<int>(g);

  if (group_element_cdf[g].empty())
  {
    new_particle.alive = false;
    return new_particle;
  }

  //======================================== Sample element
  size_t elem = std::lower_bound(
               group_element_cdf[g].begin(),
               group_element_cdf[g].end(),
               rng->Rand()) - group_element_cdf[g].begin();

  auto& src_element = group_sources[g][elem].second;

  //======================================== Sample position
  new_particle.pos = src_element.SampleRandomPosition(*rng);

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

  new_particle.cur_cell_local_id  = src_element.ParentCellLocalID();
  new_particle.cur_cell_global_id = src_element.ParentCellGlobalID();

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::Solver::RayTraceMethod::UNCOLLIDED;

  return new_particle;
}

////###################################################################
///**Gets the relative source strength accross all processors.*/
//double mcpartra::MaterialSource::GetParallelRelativeSourceWeight()
//{
//  double local_total_source_weight = 0.0;
//
//  for (auto& grp_src : group_sources)
//    for (auto& src_element : grp_src)
//      local_total_source_weight += src_element.first;
//
//  double global_total_source_weight = 0.0;
//  MPI_Allreduce(&local_total_source_weight,  //sendbuf
//                &global_total_source_weight, //recvbuf
//                1,                           //recvcount
//                MPI_DOUBLE,                  //datatype
//                MPI_SUM,                     //operation
//                MPI_COMM_WORLD);             //communicator
//
//  return local_total_source_weight/global_total_source_weight;
//}