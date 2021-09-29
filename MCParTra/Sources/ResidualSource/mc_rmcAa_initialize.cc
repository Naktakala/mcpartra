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

  size_t total_num_faces = 0;
  for (const auto& cell : grid->local_cells)
    for (const auto& face : cell.faces)
      ++total_num_faces;

  //============================================= Initialize data vectors
  //                                              (for efficiency)
  R_abs_cellk_interior.assign(num_local_cells * num_groups, 0.0);
  residual_info_cell_interiors.resize(num_groups, VecRCellInterior(num_local_cells));
  residual_info_cell_bndry_faces.resize(num_groups);
  for (auto& vec_R_cell_face : residual_info_cell_bndry_faces)
    vec_R_cell_face.reserve(total_num_faces);

  //============================================= Sample each cell
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
    size_t num_nodes = pwl->GetCellNumNodes(cell);
    const double V = cell_FV_view->volume;

    for (size_t g=0; g < num_groups; ++g)
    {
      MaterialData mat_data;
      PopulateMaterialData(cell.material_id, g, mat_data);

      auto siga = mat_data.siga;
      auto Q    = mat_data.Q;

      double cell_average_rstar_absolute=0.0;
      double cell_maximum_rstar_absolute=-1.0e-32;
      int num_points = 1000; //Number of points to sample
      for (int i=0; i < num_points; ++i)
      {
        auto x_i   = GetRandomPositionInCell(rng, cell_geom_info[k]);
        auto omega = SampleRandomDirection(rng);

        cell_pwl_view->ShapeValues(x_i, shape_values);
        cell_pwl_view->GradShapeValues(x_i,grad_shape_values);

        double phi = GetResidualFFPhi(shape_values, num_nodes, k, g);

        auto grad_phi = GetResidualFFGradPhi(
          grad_shape_values, num_nodes, cell.local_id, g);

        double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi));

        //==================================== Contribute to avg
        cell_average_rstar_absolute += std::fabs(r);
        cell_maximum_rstar_absolute = std::fmax(cell_maximum_rstar_absolute,
                                                std::fabs(r));
      }//for i
      cell_average_rstar_absolute /= num_points;

      R_abs_cellk_interior[k] = cell_average_rstar_absolute * FOUR_PI * V;

      {
        RCellInterior& rcellinterior =
          residual_info_cell_interiors[g][cell.local_id];
        rcellinterior.cell_local_id          = cell.local_id;
        rcellinterior.maximum_rstar_absolute = cell_maximum_rstar_absolute;
        rcellinterior.Rstar_absolute         = R_abs_cellk_interior[k];
      }

      {
        auto rcellinterior = std::make_unique<RCellInterior>();
        rcellinterior->cell_local_id          = cell.local_id;
        rcellinterior->maximum_rstar_absolute = cell_maximum_rstar_absolute;
        rcellinterior->Rstar_absolute         = R_abs_cellk_interior[k];

        residual_info_elements.push_back(std::move(rcellinterior));
      }

    }//for g
  }//for cell

  //============================================= Sample surfaces
  chi_log.Log(LOG_0) << "Integrating surface source.";
  for (auto& cell : grid->local_cells)
  {
    const uint64_t k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);
    size_t num_nodes = pwl->GetCellNumNodes(cell);

    unsigned int f=0;
    for (auto& face : cell.faces)
    {
      double A_f = cell_FV_view->face_area[f];
      auto&  n   = face.normal;

      for (size_t g=0; g < num_groups; ++g)
      {
        double face_average_rstar_absolute=0.0;
        double face_maximum_rstar_absolute = -1.0e32;
        int num_points = (face.has_neighbor)? 0 : 1000;
        for (int i=0; i<num_points; ++i)
        {
          auto x_i = GetRandomPositionOnCellFace(rng, cell_geom_info[k], f);

          cell_pwl_view->ShapeValues(x_i, shape_values);

          double phi = GetResidualFFPhi(shape_values, num_nodes, k, g);

          double phi_N = phi;
          if (not face.has_neighbor)
            phi_N = 0.0; //TODO: Specialize for bndries

          double r = (1.0/FOUR_PI)*(phi_N - phi);

          //==================================== Contribute to avg
          face_average_rstar_absolute += std::fabs(r);
          face_maximum_rstar_absolute  = std::fmax(face_maximum_rstar_absolute, std::fabs(r));
        }//for i
        face_average_rstar_absolute /= std::max(1, num_points);

        {
          RCellFace rcellface;
          rcellface.cell_local_id          = cell.local_id;
          rcellface.ass_face               = f;
          rcellface.maximum_rstar_absolute = face_maximum_rstar_absolute;
          rcellface.Rstar_absolute         = face_average_rstar_absolute * A_f * M_PI;

          residual_info_cell_bndry_faces[g].push_back(rcellface);
        }

        {
          auto rcellface = std::make_unique<RCellFace2>();

          rcellface->cell_local_id          = cell.local_id;
          rcellface->ass_face               = f;
          rcellface->maximum_rstar_absolute = face_maximum_rstar_absolute;
          rcellface->Rstar_absolute         = face_average_rstar_absolute * A_f * M_PI;

          residual_info_elements.push_back(std::move(rcellface));
        }
      }//for g
      ++f;
    }//for face
  }//for cell
  chi_log.Log(LOG_0) << "Preparing CDFs.";

  //============================================= Integrate sources locally
  R_abs_domain_interior.assign(num_groups, 0.0);
  R_abs_domain_surface.assign(num_groups, 0.0);
  R_abs_domain_total.assign(num_groups, 0.0);
  for (size_t g=0; g < num_groups; ++g)
  {
    R_abs_domain_interior[g] = 0.0;
    for (const auto& cell_r_info : residual_info_cell_interiors[g])
      R_abs_domain_interior[g] += cell_r_info.Rstar_absolute;

    chi_log.Log(LOG_ALL) << "Total local interior source: group " << g << " "
                         << R_abs_domain_interior[g];


    R_abs_domain_surface[g] = 0.0;
    for (auto& rcellface : residual_info_cell_bndry_faces[g])
      R_abs_domain_surface[g] += rcellface.Rstar_absolute;

    chi_log.Log(LOG_ALL) << "Total local surface source: group " << g << " "
                         << R_abs_domain_surface[g];

    R_abs_domain_total[g] = R_abs_domain_interior[g] +
                            R_abs_domain_surface[g];
  }//for g

  //============================================= Integrate residual sources
  //                                              globally
  std::vector<double> R_abs_globaldomain_interior(num_groups, 0.0);
  std::vector<double> R_abs_globaldomain_surface(num_groups, 0.0);

  double sumG_R_abs_localdomain = 0.0;
  double sumG_R_abs_globaldomain = 0.0;

  for (size_t g=0; g < num_groups; ++g)
  {
    MPI_Allreduce(&R_abs_domain_interior[g],
                  &R_abs_globaldomain_interior[g],
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);

    MPI_Allreduce(&R_abs_domain_surface[g],
                  &R_abs_globaldomain_surface[g],
                  1,
                  MPI_DOUBLE,
                  MPI_SUM,
                  MPI_COMM_WORLD);

    chi_log.Log(LOG_0) << "Total interior source: " << R_abs_globaldomain_interior[g];
    chi_log.Log(LOG_0) << "Total surface source: " << R_abs_globaldomain_surface[g];

    sumG_R_abs_localdomain += R_abs_domain_interior[g] +
                              R_abs_domain_surface[g];
    sumG_R_abs_globaldomain += R_abs_globaldomain_interior[g] +
                               R_abs_globaldomain_surface[g];
  }//for g

  SetSourceRates(sumG_R_abs_localdomain, sumG_R_abs_globaldomain);

  //============================================= Compute interior local cdf
  R_abs_domain_interior_cdf.resize(num_groups);
  for (size_t g=0; g < num_groups; ++g)
  {
    R_abs_domain_interior_cdf[g].assign(num_local_cells, 0.0);
    double running_total = 0.0;
    for (uint64_t c = 0; c < num_local_cells; ++c)
    {
      running_total += residual_info_cell_interiors[g][c].Rstar_absolute;
      R_abs_domain_interior_cdf[g][c] = running_total / R_abs_domain_interior[g];
    }
  }// for g

  //============================================= Compute surface local cdf
  R_abs_domain_surface_cdf.resize(num_groups);
  for (size_t g=0; g < num_groups; ++g)
  {
    R_abs_domain_surface_cdf[g].resize(residual_info_cell_bndry_faces[g].size(), 0.0);
    double running_total = 0.0;
    for (size_t cf = 0; cf < residual_info_cell_bndry_faces[g].size(); ++cf)
    {
      auto& rcellface = residual_info_cell_bndry_faces[g][cf];
      running_total += rcellface.Rstar_absolute;
      R_abs_domain_surface_cdf[g][cf] = running_total /
                                        R_abs_domain_surface[g];
    }
  }//for g

  //============================================= Compute groupwise
  //                                              interior-surface cdf
  interior_vs_surface_cdf.resize(num_groups, std::vector<double>(2, 0.0));
  {
    for (size_t g=0; g<num_groups; ++g)
    {
      interior_vs_surface_cdf[g][0] = R_abs_domain_interior[g] /
                                      (R_abs_domain_interior[g] + R_abs_domain_surface[g]);
      interior_vs_surface_cdf[g][1] = 1.0;
    }
  }

  //============================================= Compute group cdf
  group_cdf.assign(num_groups, 0.0);
  {
    double running_total = 0.0;
    for (size_t g=0; g<num_groups; ++g)
    {
      running_total += R_abs_domain_interior[g] +
                       R_abs_domain_surface[g];
      group_cdf[g] = running_total / sumG_R_abs_localdomain;
    }
  }
  chi_log.Log(LOG_0) << "Done initializing Residual Sources";

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