#include "mc_bndry_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Cell/cell.h"

#include "chi_log.h"

extern ChiLog& chi_log;


//###################################################################
/**Initializes references to necessary mesh elements.
 *
 * Since this is a boundary source the major difficulty here is
 * to get the sampling direction right. Sampled directions need
 * to be rotation to the opposite of the outward pointing boundary
 * normals. For this we assemble a rotation matrix for each face
 * that needs to be sampled.*/
void mcpartra::BoundaryIsotropicSource::
  Initialize(chi_mesh::MeshContinuumPtr&    ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
             const std::vector<CellGeometryData>& ref_cell_geometry_info)
{
  chi_log.Log(LOG_0) << "Initializing Boundary Source";

  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;

  const int ALL_BOUNDRIES = -1;

  //============================================= Make checks
  if (mg_isotropic_strength.size() != num_groups)
    throw std::invalid_argument("Multigroup isotropic source strength supplied "
                                "to isotropic boundary source does not have "
                                "a compatible number of groups.");

  //============================================= Initialize mg containers
  group_elements.clear();
  group_elements.resize(num_groups);

  group_element_cdf.clear();
  group_element_cdf.resize(num_groups);

  group_cdf.clear();
  group_cdf.assign(num_groups, 0.0);

  //============================================= Build surface src patchess
  // The surface points
  double IntA_Q_total_local = 0.0;
  std::vector<double> IntA_Q_g(num_groups, 0.0);
  for (auto& cell : grid->local_cells)
  {
    auto fv_view = fv_sdm->MapFeView(cell.local_id);

    int f=0;
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor) {++f; continue;}

      //Determine if face will be sampled
      bool face_is_a_source = false;
      if (ref_bndry == ALL_BOUNDRIES) face_is_a_source = true;
      if ((ref_bndry != ALL_BOUNDRIES) and (ref_bndry == face.neighbor_id))
        face_is_a_source = true;

      if (face_is_a_source)
      {
        double face_area = fv_view->face_area[f];

        for (size_t g=0; g<num_groups; ++g)
        {
          double source_strength = face_area * M_PI * mg_isotropic_strength[g];
          if (source_strength < 1.0e-12) continue;

          auto face_elements =
            mcpartra::GetCellSurfaceSourceElements(cell, face, grid);

          for (auto& element : face_elements)
          {
            double element_strength = source_strength*element.Area()/face_area;
            if (element_strength < 1.0e-12) continue;
            group_elements[g].emplace_back(element_strength, element);
            IntA_Q_g[g] += element_strength;
          }//for element
        }//for g

      }
      ++f;
    }//for face
  }//for cell

  for (auto val : IntA_Q_g)
      IntA_Q_total_local += val;

  //======================================== Compute normalized pdf used during
  //                                         biasing

  //============================================= Set source strength for
  //                                              solver
  {
    double IntA_Q_total_globl = 0.0;
    MPI_Allreduce(&IntA_Q_total_local,     // sendbuf
                  &IntA_Q_total_globl,     // recvbuf
                  1, MPI_DOUBLE,           // count + datatype
                  MPI_SUM,                 // operation
                  MPI_COMM_WORLD);         // communicator

    SetSourceRates(IntA_Q_total_local,IntA_Q_total_globl);
  }

  //============================================= Build element-wise CDFs
  for (size_t g=0; g<num_groups; ++g)
  {
    const size_t num_elements = group_elements[g].size();
    group_element_cdf[g].assign(num_elements, 0.0);

    double cumulative_value = 0.0;
    size_t ell = 0;
    for (const auto& src_element : group_elements[g])
    {
      cumulative_value += src_element.first;
      group_element_cdf[g][ell] = cumulative_value/IntA_Q_g[g];
      ++ell;
    }
  }//for g

  //============================================= Build group-wise CDFs
  {
    double cumulative_value = 0.0;
    for (size_t g=0; g<num_groups; ++g)
    {
      cumulative_value += IntA_Q_g[g];
      group_cdf[g] = cumulative_value/IntA_Q_total_local;
    }//for g
  }

}

//###################################################################
/**Source routine for boundary source.*/
mcpartra::Particle mcpartra::BoundaryIsotropicSource::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  const std::string fname = __FUNCTION__;

  mcpartra::Particle new_particle;

  //======================================== Sample group
  size_t g = SampleCDF(group_cdf, rng);

  new_particle.egrp = static_cast<int>(g);

  if (group_element_cdf[g].empty())
  {new_particle.alive = false; return new_particle;}

  //======================================== Sample element
  size_t elem = SampleCDF(group_element_cdf[g], rng);

  //======================================== Get references
  const auto& src_element = group_elements[g].at(elem);
  const auto& element = src_element.second;

  //======================================== Sample position
  new_particle.pos = element.SampleRandomPosition(rng);

  //======================================== Sample direction
  new_particle.dir = mcpartra::RandomCosineLawDirection(rng, -1*element.Normal());

  //======================================== Determine weight
  new_particle.w = 1.0;

  new_particle.cur_cell_global_id = element.ParentCellGlobalID();
  new_particle.cur_cell_local_id  = element.ParentCellLocalID();

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::RayTraceMethod::UNCOLLIDED;

  return new_particle;
}