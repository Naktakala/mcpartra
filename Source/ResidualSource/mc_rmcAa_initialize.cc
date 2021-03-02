#include "mc_rmcA_source.h"

#include "../../Solver/solver_montecarlon.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualSourceA::
ResidualSourceA(std::shared_ptr<chi_physics::FieldFunction> in_resid_ff,
                bool use_uniform_sampling) :
  sample_uniformly(use_uniform_sampling)
{
  type_index = SourceTypes::RESIDUAL_TYPE_A;
  resid_ff = in_resid_ff;
}

//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void chi_montecarlon::ResidualSourceA::
Initialize(chi_mesh::MeshContinuumPtr ref_grid,
           std::shared_ptr<SpatialDiscretization_FV> ref_fv_sdm,
           chi_montecarlon::Solver* ref_solver)
{
  chi_log.Log(LOG_0) << "Initializing Residual3 Sources";

  const double FOUR_PI   = 4.0*M_PI;
  grid                   = ref_grid;
  fv_sdm                 = ref_fv_sdm;
  this->ref_solver       = ref_solver;
  auto& rng              = ref_solver->rng0;
  size_t num_local_cells = ref_solver->grid->local_cells.size();
  auto& pwl              = ref_solver->pwl;
  auto& fv               = ref_solver->fv;

  //============================================= Build cell composition data
  //This info is basically constituent side triangles/tetrahedrons
  BuildCellVolInfo(ref_grid,ref_fv_sdm);

  //============================================= Transform tilde_phi to
  //                                              tilde_phi_star
  RemoveFFDiscontinuities();

  //============================================= Initialize data vectors
  //                                              (for efficiency)
  r_abs_cellk_interior_average.      clear();
  r_cellk_interior_max.              clear();
  r_cellk_interior_min.              clear();
  R_abs_cellk_interior.              clear();
//  R_abs_cellk_surface.               clear();

  r_abs_cellk_interior_average.      resize(num_local_cells, 0.0);
  r_cellk_interior_max.              resize(num_local_cells, -1.0e32);
  r_cellk_interior_min.              resize(num_local_cells,  1.0e32);
  R_abs_cellk_interior.              resize(num_local_cells, 0.0);
//  R_abs_cellk_surface.               resize(num_local_cells, 0.0);

  //============================================= Sample each cell
  chi_log.Log(LOG_0) << "Integrating cell source.";

  //Predefine N_i and grad_N_i vectors to prevent
  //unnecessary memory reallocation.
  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;
  MaterialData                   mat_data;

  for (auto& cell : ref_solver->grid->local_cells)
  {
    const int k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);

    int group_g = 0;

    //====================================== Get material properties
    PopulateMaterialData(cell.material_id,group_g,mat_data);

    auto siga = mat_data.siga;
    auto Q    = mat_data.Q;

    //====================================== Determine absolute Residual Interior
    int N_p = 1000; //Number of points to sample
    for (int i=0; i < N_p; ++i)
    {
      auto x_i = GetRandomPositionInCell(rng, cell_geometry_info[k]);

      chi_mesh::Vector3 omega = RandomDirection(rng);

      cell_pwl_view->ShapeValues(x_i, shape_values);
      cell_pwl_view->GradShapeValues(x_i,grad_shape_values);

      double phi = GetResidualFFPhi(shape_values, cell_pwl_view->num_nodes, k, group_g);

      auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
                                           cell_pwl_view->num_nodes,
                                           cell.local_id,
                                           0);

      double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi));

      //==================================== Contribute to avg
      r_abs_cellk_interior_average[k] += std::fabs(r);
      r_cellk_interior_max[k] = std::fmax(r_cellk_interior_max[k], std::fabs(r));
      r_cellk_interior_min[k] = std::fmin(r_cellk_interior_min[k], r);
    }//for i
    r_abs_cellk_interior_average[k]  /= N_p;

    R_abs_cellk_interior[k] = FOUR_PI * cell_FV_view->volume *
                              r_abs_cellk_interior_average[k];

    domain_volume     += cell_FV_view->volume;

    if (std::fabs(r_abs_cellk_interior_average[k]) > 1.0e-16)
      source_volume += cell_FV_view->volume;
  }//for cell

  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Sample surfaces
  chi_log.Log(LOG_0) << "Integrating surface source.";
  for (auto& cell : ref_solver->grid->local_cells)
  {
    const int k = cell.local_id;
    auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);
    auto cell_FV_view  = fv->MapFeView(cell.local_id);

    int group_g = 0;

    int f=-1;
    for (auto& face : cell.faces)
    {
      ++f;
      double A_f = cell_FV_view->face_area[f];
      auto&  n   = face.normal;

      RCellFace rcellface;
      rcellface.cell_local_id = cell.local_id;
      rcellface.ass_face      = f;
      rcellface.area          = A_f;

      double face_ave_surface_pstar=0.0;
      int num_points = (face.has_neighbor)? 0 : 1000;
      for (int i=0; i<num_points; ++i)
      {
        auto x_i = GetRandomPositionOnCellSurface(rng, cell_geometry_info[k], f);

        cell_pwl_view->ShapeValues(x_i, shape_values);

        double phi = GetResidualFFPhi(shape_values, cell_pwl_view->num_nodes, k, group_g);

        double phi_N = phi;
        if (not face.has_neighbor)
          phi_N = 0.0; //TODO: Specialize for bndries

        double r = (1.0/FOUR_PI)*(phi_N - phi);

        //==================================== Contribute to avg
        face_ave_surface_pstar += std::fabs(r);
        rcellface.maximum = std::fmax(rcellface.maximum,std::fabs(r));
        rcellface.minimum = std::fmin(rcellface.minimum,r);
      }//for i
      face_ave_surface_pstar /= std::max(1,num_points);

      rcellface.average_rstar = face_ave_surface_pstar;

      r_abs_cellk_facef_surface_average.push_back(rcellface);
    }//for face
  }//for cell

  //============================================= Integrate sources locally
  R_abs_localdomain_interior = 0.0;
  for (double v : R_abs_cellk_interior)
    R_abs_localdomain_interior += v;

  chi_log.Log(LOG_ALL) << "Total local interior source: " << R_abs_localdomain_interior;

  R_abs_localdomain_surface = 0.0;
  for (auto& rcellface : r_abs_cellk_facef_surface_average)
    R_abs_localdomain_surface += rcellface.average_rstar * rcellface.area * M_PI;

  chi_log.Log(LOG_ALL) << "Total local surface source: " << R_abs_localdomain_surface;

  //============================================= Integrate residual sources
  //                                              globally
  MPI_Allreduce(&R_abs_localdomain_interior,
                &R_abs_globaldomain_interior,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

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
    for (int c = 0; c < num_local_cells; ++c)
    {
      running_total += R_abs_cellk_interior[c];
      domain_cdf[c] = running_total / R_abs_localdomain_interior;
    }
  }

  //============================================= Compute surface local cdf
  {
    surface_cdf.clear();
    surface_cdf.resize(r_abs_cellk_facef_surface_average.size(), 0.0);
    double running_total = 0.0;
    for (int cf = 0; cf < r_abs_cellk_facef_surface_average.size(); ++cf)
    {
      auto& rcellface = r_abs_cellk_facef_surface_average[cf];
      running_total += rcellface.average_rstar*rcellface.area*M_PI;
      surface_cdf[cf] = running_total / R_abs_localdomain_surface;
    }
  }
  chi_log.Log(LOG_0) << "Done initializing Residual Sources";

  ref_solver->source_normalization = R_abs_globaldomain_interior +
                                     R_abs_globaldomain_surface;

}