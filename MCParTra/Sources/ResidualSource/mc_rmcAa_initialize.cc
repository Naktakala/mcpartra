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

  //============================================= Initialize data vectors
  //                                              (for efficiency)
  R_abs_cellk_interior.assign(num_local_cells, 0.0);

  //============================================= Sample each cell
  chi_log.Log(LOG_0) << "Integrating cell source.";

  //Predefine N_i and grad_N_i vectors to prevent
  //unnecessary memory reallocation.
  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  residual_info_cell_interiors.resize(num_local_cells);
  for (const auto& cell : grid->local_cells)
  {
    const uint64_t k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);
    size_t num_nodes = pwl->GetCellNumNodes(cell);
    const double V = cell_FV_view->volume;

    int group_g = 0;

    //====================================== Get material properties
    MaterialData mat_data;
    PopulateMaterialData(cell.material_id,group_g,mat_data);

    auto siga = mat_data.siga;
    auto Q    = mat_data.Q;

    //====================================== Determine absolute Residual Interior
    double cell_average_rstar_absolute=0.0;
    double cell_maximum_rstar_absolute=-1.0e-32;
    int num_points = 1000; //Number of points to sample
    for (int i=0; i < num_points; ++i)
    {
      auto x_i   = GetRandomPositionInCell(rng, cell_geom_info[k]);
      auto omega = SampleRandomDirection(rng);

      cell_pwl_view->ShapeValues(x_i, shape_values);
      cell_pwl_view->GradShapeValues(x_i,grad_shape_values);

      double phi = GetResidualFFPhi(shape_values, num_nodes, k, group_g);

      auto grad_phi = GetResidualFFGradPhi(
        grad_shape_values, num_nodes, cell.local_id, group_g);

      double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi));

      //==================================== Contribute to avg
      cell_average_rstar_absolute += std::fabs(r);
      cell_maximum_rstar_absolute = std::fmax(cell_maximum_rstar_absolute,
                                              std::fabs(r));
    }//for i
    cell_average_rstar_absolute /= num_points;

    R_abs_cellk_interior[k] = cell_average_rstar_absolute * FOUR_PI * V;

    RCellInterior& rcellinterior = residual_info_cell_interiors[cell.local_id];
    rcellinterior.cell_local_id          = cell.local_id;
    rcellinterior.average_rstar_absolute = cell_average_rstar_absolute;
    rcellinterior.maximum_rstar_absolute = cell_maximum_rstar_absolute;
    rcellinterior.Rstar_absolute         = R_abs_cellk_interior[k];
  }//for cell

  //============================================= Sample surfaces
  chi_log.Log(LOG_0) << "Integrating surface source.";
  for (auto& cell : grid->local_cells)
  {
    const uint64_t k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);

    int group_g = 0;

    int f=-1;
    for (auto& face : cell.faces)
    {
      ++f;
      double A_f = cell_FV_view->face_area[f];
      auto&  n   = face.normal;


      double face_average_rstar_absolute=0.0;
      double face_maximum_rstar_absolute = -1.0e32;
      int num_points = (face.has_neighbor)? 0 : 1000;
      for (int i=0; i<num_points; ++i)
      {
        auto x_i = GetRandomPositionOnCellFace(rng, cell_geom_info[k], f);

        cell_pwl_view->ShapeValues(x_i, shape_values);

        double phi = GetResidualFFPhi(shape_values, cell_pwl_view->num_nodes, k, group_g);

        double phi_N = phi;
        if (not face.has_neighbor)
          phi_N = 0.0; //TODO: Specialize for bndries

        double r = (1.0/FOUR_PI)*(phi_N - phi);

        //==================================== Contribute to avg
        face_average_rstar_absolute += std::fabs(r);
        face_maximum_rstar_absolute  = std::fmax(face_maximum_rstar_absolute, std::fabs(r));

      }//for i
      face_average_rstar_absolute /= std::max(1, num_points);

      RCellFace rcellface;
      rcellface.cell_local_id          = cell.local_id;
      rcellface.ass_face               = f;
      rcellface.maximum_rstar_absolute = face_maximum_rstar_absolute;
      rcellface.Rstar_absolute         = face_average_rstar_absolute * A_f * M_PI;

      residual_info_cell_bndry_faces.push_back(rcellface);
    }//for face
  }//for cell

  //============================================= Integrate sources locally
  R_abs_localdomain_interior = 0.0;
  for (const auto& cell_r_info : residual_info_cell_interiors)
    R_abs_localdomain_interior += cell_r_info.Rstar_absolute;

  chi_log.Log(LOG_ALL) << "Total local interior source: "
                       << R_abs_localdomain_interior;

  R_abs_localdomain_surface = 0.0;
  for (auto& rcellface : residual_info_cell_bndry_faces)
    R_abs_localdomain_surface += rcellface.Rstar_absolute;

  chi_log.Log(LOG_ALL) << "Total local surface source: "
                       << R_abs_localdomain_surface;

  //============================================= Integrate residual sources
  //                                              globally
  double R_abs_globaldomain_interior;
  MPI_Allreduce(&R_abs_localdomain_interior,
                &R_abs_globaldomain_interior,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  double R_abs_globaldomain_surface;
  MPI_Allreduce(&R_abs_localdomain_surface,
                &R_abs_globaldomain_surface,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi_log.Log(LOG_0) << "Total interior source: " << R_abs_globaldomain_interior;
  chi_log.Log(LOG_0) << "Total surface source: " << R_abs_globaldomain_surface;

  //============================================= Compute interior local cdf
  {
    domain_cdf.clear();
    domain_cdf.resize(num_local_cells, 0.0);
    double running_total = 0.0;
    for (uint64_t c = 0; c < num_local_cells; ++c)
    {
      running_total += R_abs_cellk_interior[c];
      domain_cdf[c] = running_total / R_abs_localdomain_interior;
    }
  }

  //============================================= Compute surface local cdf
  {
    surface_cdf.clear();
    surface_cdf.resize(residual_info_cell_bndry_faces.size(), 0.0);
    double running_total = 0.0;
    for (size_t cf = 0; cf < residual_info_cell_bndry_faces.size(); ++cf)
    {
      auto& rcellface = residual_info_cell_bndry_faces[cf];
      running_total += rcellface.Rstar_absolute;
      surface_cdf[cf] = running_total / R_abs_localdomain_surface;
    }
  }
  chi_log.Log(LOG_0) << "Done initializing Residual Sources";

  SetSourceRates(R_abs_localdomain_interior + R_abs_localdomain_surface,
                 R_abs_globaldomain_interior + R_abs_globaldomain_surface);

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