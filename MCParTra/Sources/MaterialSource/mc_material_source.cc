#include "mc_material_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Initialize material source.*/
void mcpartra::MaterialSource::
  Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
             const std::vector<CellGeometryData>& ref_cell_geometry_info)
{
  chi_log.Log(LOG_0) << "Initializing Material Sources";

  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;

  size_t num_materials = chi_physics_handler.material_stack.size();

  //============================================= List cells with mat-sources
  std::map<uint64_t, chi_mesh::Cell*> mat_src_cells;
  matid_has_q_flags.assign(num_materials, false);
  for (auto& cell : ref_grid->local_cells)
  {
    auto mat_id = static_cast<int>(cell.material_id);
    auto mat = chi_physics_handler.material_stack[mat_id];
    for (const auto& prop : mat->properties)
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        mat_src_cells.insert(std::make_pair(cell.local_id,&cell));
        cell_elements.insert(
          std::make_pair(cell.local_id, GetCellVolumeSourceElements(cell,grid)));

        matid_has_q_flags[mat_id] = true;
        matid_q_map[mat_id] = std::static_pointer_cast<IsoMGSrc>(prop);
        break;
      }
  }

  //============================================= Determine group-wise
  //                                              source weights
  std::vector<double> IntV_Q_g(ref_solver.num_groups, 0.0);
  group_sources.resize(ref_solver.num_groups);
  for (auto& cell_index_ptr_pair : mat_src_cells)
  {
    auto& cell = *cell_index_ptr_pair.second;
    auto fv_view = ref_fv_sdm->MapFeView(cell.local_id);

    auto mat_id = static_cast<int>(cell.material_id);
    auto mat = chi_physics_handler.material_stack[mat_id];

    if (matid_has_q_flags[mat_id])
    {
      auto src = matid_q_map[mat_id];
      for (size_t g=0; g<src->source_value_g.size(); ++g)
      {
        double Q_g = src->source_value_g[g];

        if (not (std::fabs(Q_g) > 0.0)) continue;

        IntV_Q_g[g] += fv_view->volume*Q_g;

        auto& v_elements = cell_elements.at(cell.local_id);
        for (auto& element : v_elements)
          group_sources[g].emplace_back(Q_g * element.Volume(), element);
      }//for g
    }//if has source
  }//for cell

  //============================================= Checking group sources integrity
  for (const auto& group_source : group_sources)
    for (const auto& source_item : group_source)
      if (source_item.second.TypeIndex() == 0)
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
  group_cdf.resize(ref_solver.num_groups, 0.0);
  double running_total = 0.0;
  for (size_t g=0; g<ref_solver.num_groups; ++g)
  {
    running_total += IntV_Q_g[g];
    group_cdf[g] = (IntV_Q_total_local > 0.0)?
                    running_total / IntV_Q_total_local : 0.0;
  }

  //============================================= Within each group, construct
  //                                              source element cdf
  group_element_cdf.clear();
  group_element_cdf.resize(ref_solver.num_groups);
  size_t g=0;
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
  const std::string fname{__FUNCTION__};

  Particle new_particle;

  //======================================== Sample group
  size_t g = SampleCDF(group_cdf, rng);

  if (g >= group_cdf.size())
    throw std::logic_error(fname + ": Error sampling group cdf.");

  new_particle.egrp = static_cast<int>(g);

  if (group_element_cdf[g].empty())
  {
    new_particle.alive = false;
    return new_particle;
  }

  //======================================== Sample element
  auto elem = static_cast<size_t>(
    std::lower_bound(
               group_element_cdf[g].begin(),
               group_element_cdf[g].end(),
               rng.Rand()) - group_element_cdf[g].begin());

  if (elem >= group_sources[g].size())
    throw std::logic_error(std::string(__FUNCTION__) +
                           ": Error sampling group_sources[" +
                           std::to_string(g) + "]");

  auto& src_element = group_sources[g][elem].second;

  //======================================== Sample position
  new_particle.pos = src_element.SampleRandomPosition(rng);

  //======================================== Sample direction
  new_particle.dir = SampleRandomDirection(rng);

  //======================================== Determine weight
  new_particle.w = 1.0;

  new_particle.cur_cell_local_id  = src_element.ParentCellLocalID();
  new_particle.cur_cell_global_id = src_element.ParentCellGlobalID();

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::RayTraceMethod::UNCOLLIDED;

  //======================================== Determine moment indices
  const auto phi_theta = OmegaToPhiThetaSafe(new_particle.dir);
  const double& phi   = phi_theta.first;
  const double& theta = phi_theta.second;

  size_t num_moments = ref_solver.num_moments;
  for (size_t m=1; m<num_moments; ++m)
  {
    const auto& ell_em = ref_solver.m_to_ell_em_map[m];
    auto ell = static_cast<unsigned int>(ell_em.first);
    int em           = ell_em.second;

    new_particle.moment_values[m] = chi_math::Ylm(ell,em,phi,theta);
  }

  return new_particle;
}