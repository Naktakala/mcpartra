#include "mc_material_source.h"

#include "../../Solver/solver_montecarlon.h"

#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include "ChiMesh/Cell/cell_polyhedron.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

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

          auto& v_elements = cell_elements.at(cell.local_id);
          for (auto& element : v_elements)
            group_sources[g].emplace_back(Q_g * element.Volume(), &element);
        }//for g
      }//if src prop
  }//for cell

  //============================================= Checking group sources integrity
  for (const auto& group_source : group_sources)
    for (const auto& source_item : group_source)
      if (source_item.second->TypeIndex() == 0)
        throw std::logic_error(std::string(__FUNCTION__) +
                               ": Group source integrity failed.");

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
    group_cdf[g] = (IntV_Q_total_local > 0.0)?
                    running_total / IntV_Q_total_local : 0.0;
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
      group_element_cdf[g][elem] = (IntV_Q_g[g] > 0.0)?
                                   elem_running_total/IntV_Q_g[g] : 0.0;
      ++elem;
    }
    ++g;
  }

}


//###################################################################
/**Creates a source particle.*/
mcpartra::Particle mcpartra::MaterialSource::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  Particle new_particle;

  //======================================== Sample group
  size_t g = std::lower_bound(
            group_cdf.begin(),
            group_cdf.end(),
            rng.Rand()) - group_cdf.begin();

  if (g >= group_cdf.size())
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Error sampling group cdf.");

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
               rng.Rand()) - group_element_cdf[g].begin();

  if (elem >= group_sources[g].size())
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Error sampling group_sources[" +
                           std::to_string(g) + "]");

  auto& src_element = group_sources[g][elem].second;

  //======================================== Sample position
  new_particle.pos = src_element->SampleRandomPosition(rng);

  //======================================== Sample direction
  new_particle.dir = SampleRandomDirection(rng);

  //======================================== Determine weight
  new_particle.w = 1.0;

  new_particle.cur_cell_local_id  = src_element->ParentCellLocalID();
  new_particle.cur_cell_global_id = src_element->ParentCellGlobalID();

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::Solver::RayTraceMethod::UNCOLLIDED;

  //======================================== Determine moment indices
  int num_moments = ref_solver.num_moms;
  for (int m=1; m<num_moments; ++m)
  {
    const auto& ell_em = ref_solver.m_to_ell_em_map[m];
    int ell = ell_em.first;
    int em = ell_em.second;

    auto phi_theta = OmegaToPhiThetaSafe(new_particle.dir);
    double& phi   = phi_theta.first;
    double& theta = phi_theta.second;

    new_particle.moment_values[m] = chi_math::Ylm(ell,em,phi,theta);
  }

  return new_particle;
}