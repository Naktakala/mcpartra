#include "mc_rmcA_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;


//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void mcpartra::ResidualSourceA::
  Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
             const std::vector<CellGeometryData>& ref_cell_geometry_info)
{
  typedef chi_mesh::Vector3 Vec3;
  chi_log.Log(LOG_0) << "Initializing Residual3 Sources";

  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;

  cell_geometry_info =
    std::make_unique<std::vector<CellGeometryData>>(ref_cell_geometry_info);
  auto& cell_geom_info = *cell_geometry_info;

  auto& rng              = ref_solver.rng0;
  auto& pwl              = ref_solver.pwl;
  auto& fv               = fv_sdm;

  const double FOUR_PI   = 4.0*M_PI;
  const size_t num_local_cells = grid->local_cells.size();

  //============================================= Transform tilde_phi to
  //                                              tilde_phi_star
  RemoveFFDiscontinuities();

  size_t total_num_faces = 0;
  for (const auto& cell : grid->local_cells)
    for (const auto& face : cell.faces)
      ++total_num_faces;

  //============================================= Initialize data vectors
  //                                              (for efficiency)
  group_sources.resize(num_groups);
  R_abs_cellk_interior.assign(num_local_cells * num_groups, 0.0);

  //============================================= Determine group-wise
  //                                              interior source strengths
  chi_log.Log(LOG_0) << "Integrating cell source.";

  //Predefine N_i and grad_N_i vectors to prevent
  //unnecessary memory reallocation.
  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  for (const auto& cell : grid->local_cells)
  {
    const uint64_t k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);
    const size_t num_nodes = pwl->GetCellNumNodes(cell);
    const double V = cell_FV_view->volume;

    for (size_t g=0; g < num_groups; ++g)
    {
      MaterialData mat_data;
      PopulateMaterialData(cell.material_id, g, mat_data);

      auto siga = mat_data.siga;
      auto Q    = mat_data.Q;

      const VecDbl& nodal_phi = GetResidualFFPhiAtNodes(cell, num_nodes, 0, g);

      double cell_average_rstar_absolute=0.0;
      double cell_maximum_rstar_absolute=-1.0e-32;
      int num_points = 1000; //Number of points to sample
      for (int i=0; i < num_points; ++i)
      {
        auto x_i   = GetRandomPositionInCell(rng, cell_geom_info[k]);
        auto omega = SampleRandomDirection(rng);

        cell_pwl_view->ShapeValues(x_i, shape_values);
        cell_pwl_view->GradShapeValues(x_i,grad_shape_values);

        double phi = GetPhiH(shape_values, nodal_phi, num_nodes);
        Vec3   grad_phi = GetGradPhiH(grad_shape_values, nodal_phi, num_nodes);

        double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi));

        cell_average_rstar_absolute += std::fabs(r);
        cell_maximum_rstar_absolute = std::fmax(cell_maximum_rstar_absolute,
                                                std::fabs(r));
      }//for i
      cell_average_rstar_absolute /= num_points;

      R_abs_cellk_interior[k] = cell_average_rstar_absolute * FOUR_PI * V;

      auto rcellinterior = std::make_unique<RCellInterior>();
      rcellinterior->cell_local_id          = cell.local_id;
      rcellinterior->maximum_rstar_absolute = cell_maximum_rstar_absolute;
      rcellinterior->Rstar_absolute         = R_abs_cellk_interior[k];

      group_sources[g].push_back(std::move(rcellinterior));
    }//for g
  }//for cell

  //============================================= Detemine group-wise
  //                                              face source strengths
  chi_log.Log(LOG_0) << "Integrating surface source.";
  for (auto& cell : grid->local_cells)
  {
    const uint64_t k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);
    const size_t num_nodes = pwl->GetCellNumNodes(cell);

    unsigned int f=0;
    for (auto& face : cell.faces)
    {
      double A_f = cell_FV_view->face_area[f];
      auto&  n   = face.normal;

      for (size_t g=0; g < num_groups; ++g)
      {
        const VecDbl& nodal_phi = GetResidualFFPhiAtNodes(cell, num_nodes, 0, g);

        double face_average_rstar_absolute=0.0;
        double face_maximum_rstar_absolute = -1.0e32;
        int num_points = (face.has_neighbor)? 0 : 1000;
        for (int i=0; i<num_points; ++i)
        {
          auto x_i = GetRandomPositionOnCellFace(rng, cell_geom_info[k], f);

          cell_pwl_view->ShapeValues(x_i, shape_values);

          double phi = GetPhiH(shape_values, nodal_phi, num_nodes);

          double phi_N = phi;
          if (not face.has_neighbor)
            phi_N = 0.0; //TODO: Specialize for bndries

          double r = (1.0/FOUR_PI)*(phi_N - phi);

          //==================================== Contribute to avg
          face_average_rstar_absolute += std::fabs(r);
          face_maximum_rstar_absolute  = std::fmax(face_maximum_rstar_absolute, std::fabs(r));
        }//for i
        face_average_rstar_absolute /= std::max(1, num_points);

        auto rcellface = std::make_unique<RCellFace>();

        rcellface->cell_local_id          = cell.local_id;
        rcellface->ass_face               = f;
        rcellface->maximum_rstar_absolute = face_maximum_rstar_absolute;
        rcellface->Rstar_absolute         = face_average_rstar_absolute * A_f * M_PI;

        group_sources[g].push_back(std::move(rcellface));

        R_abs_cellk_interior[cell.local_id] += face_average_rstar_absolute * A_f * M_PI;
      }//for g
      ++f;
    }//for face
  }//for cell

  chi_log.Log(LOG_0) << "Preparing CDFs.";

  //============================================= Determine totals
  std::vector<double> IntV_Q_g(num_groups, 0.0);
  for (size_t g=0; g<num_groups; ++g)
    for (const auto& rinfo : group_sources[g])
      IntV_Q_g[g] += rinfo->Rstar_absolute;

  double IntV_Q_total_local = 0.0;
  for (auto val : IntV_Q_g) IntV_Q_total_local += val;

  //======================================== Compute normalized pdf used during
  //                                         biasing
  group_element_pdf.resize(num_groups);
  for (size_t g=0; g<num_groups; ++g)
  {
    size_t num_elems = group_sources[g].size();
    group_element_pdf[g].resize(num_elems, 0.0);
    for (size_t elem=0; elem<num_elems; ++elem)
      group_element_pdf[g][elem] =
        group_sources[g][elem]->Rstar_absolute / IntV_Q_g[g];
  }

  //============================================= Set source strength for
  //                                              solver
  {
    double IntV_Q_total_globl = 0.0;
    MPI_Allreduce(&IntV_Q_total_local,     // sendbuf
                  &IntV_Q_total_globl,     // recvbuf
                  1, MPI_DOUBLE,           // count + datatype
                  MPI_SUM,                 // operation
                  MPI_COMM_WORLD);         // communicator

    SetSourceRates(IntV_Q_total_local,IntV_Q_total_globl);
  }

  //============================================= Construct groupwise cdf
  group_cdf.clear();
  group_cdf.resize(num_groups, 0.0);
  {
    double running_total = 0.0;
    for (size_t g=0; g<num_groups; ++g)
    {
      running_total += IntV_Q_g[g];
      group_cdf[g] = (IntV_Q_total_local > 0.0) ?
        running_total / IntV_Q_total_local : 0.0;
    }//for g
  }

  //============================================= Within each group, construct
  //                                              source element cdf
  group_element_cdf.clear();
  group_element_cdf.resize(num_groups);
  {
    size_t g = 0;
    for (auto& group_source : group_sources)
    {
      group_element_cdf[g].resize(group_source.size(), 0.0);
      double elem_running_total = 0.0;
      size_t elem = 0;
      for (auto& src_element : group_source)
      {
        elem_running_total += src_element->Rstar_absolute;
        group_element_cdf[g][elem] = (IntV_Q_g[g] > 0.0) ?
          elem_running_total/IntV_Q_g[g] : 0.0;
        ++elem;
      }//for elem
      ++g;
    }//for g
  }

  chi_log.Log(LOG_0) << "Done initializing Residual Sources";

  BiasCDFs(false); //This just copies cdfs and corrections

  //============================================= Export interior source
  //                                              as FieldFunction
  auto fv_sd = std::dynamic_pointer_cast<SpatialDiscretization>(fv);
  auto R_ff = std::make_shared<chi_physics::FieldFunction>(
    "R_interior",                                 //Text name
    fv_sd,                                        //Spatial Discretization
    &R_abs_cellk_interior,                        //Data
    ref_solver.uk_man_fv,                         //Nodal variable structure
    0, 0);                                        //Reference variable and component

  R_ff->ExportToVTKFV("ZRout","R_interior");



}