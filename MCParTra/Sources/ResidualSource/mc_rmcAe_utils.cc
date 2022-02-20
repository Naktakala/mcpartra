#include "mc_rmcA_source.h"

#include "ChiMesh/Cell/cell.h"

#include "ChiMath/Quadratures/product_quadrature.h"
#include "ChiMath/Quadratures/LegendrePoly/legendrepoly.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "ChiMiscUtils/chi_misc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Populates a material information data structure from a mat-id.*/
void mcpartra::ResidualSourceA::
  PopulateMaterialData(int mat_id, size_t group_g, MaterialData &mat_data)
{
  auto xs = ref_solver.matid_xs_map2[mat_id];
  double siga = xs->sigma_a[group_g];
  double sigt = xs->sigma_t[group_g];

  double Q    = 0.0;
  if (ref_solver.matid_has_q_flags[mat_id])
  {
    auto q_prop = ref_solver.matid_q_map2[mat_id];
    Q = q_prop->source_value_g[group_g];
  }

  mat_data.siga = siga;
  mat_data.Q = Q;
}

//###################################################################
/**Obtains a field function interpolant of the flux.*/
std::vector<double> mcpartra::ResidualSourceA::
  GetResidualFFPhiAtNodes(const chi_mesh::Cell &cell,
                          size_t num_nodes,
                          size_t variable_id,
                          size_t component_id)
{
  typedef SpatialDiscretization_FE SDMFE;
  typedef chi_math::SpatialDiscretizationType SDType;
  const auto SDTypePWLD = SDType::PIECEWISE_LINEAR_DISCONTINUOUS;

  const auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != SDTypePWLD)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SDMFE>(sdm);

  const auto& uk_man = resid_ff->unknown_manager;

  std::vector<double> phi(num_nodes, 0.0);
  for (size_t n=0; n < num_nodes; n++) //node n
  {
    int64_t dof_map = pwl_sdm->MapDOFLocal(cell,
                                           n,
                                           uk_man,
                                           variable_id,
                                           component_id);
    phi[n] = (*resid_ff->field_vector_local)[dof_map];
  }//for node n

  return phi;
}


//###################################################################
/**Obtains a field function interpolant of the flux.*/
double mcpartra::ResidualSourceA::
  GetPhiH(const std::vector<double>& shape_values,
          const std::vector<double>& phi,
          size_t num_nodes)
{
  double phi_h = 0.0;
  for (size_t j=0; j < num_nodes; ++j)
    phi_h += phi[j] * shape_values[j];

  return phi_h;
}


//###################################################################
/**Obtains a field function interpolant of the angular flux.*/
double mcpartra::ResidualSourceA::
  GetPsiH(const std::vector<double>& shape_values,
          const std::vector<VecDbl>& phi,
          const chi_mesh::Vector3 omega,
          size_t num_nodes,
          const std::vector<EllEmIndices>& m_to_ell_em_map)
{
  const size_t num_moms = phi.size();

  const auto phi_theta   = OmegaToPhiThetaSafe(omega);
  const double& varphi   = phi_theta.first;
  const double& theta    = phi_theta.second;

  double psi_h = 0.0;
  for (size_t m=0; m<num_moms; ++m)
  {
    const EllEmIndices& ell_em = m_to_ell_em_map[m];
    const unsigned int ell = ell_em.ell;
    const int           em = ell_em.m;

    const double Y_ell_em = chi_math::Ylm(ell, em, varphi, theta);
    const double factor = ((2.0*ell + 1.0)/4.0/M_PI);

    for (size_t j=0; j < num_nodes; ++j)
      psi_h += factor * Y_ell_em * phi[m][j] * shape_values[j];
  }

  return psi_h;
}


//###################################################################
/**Obtains a field function interpolant of the flux-gradient.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  GetGradPhiH(
    const std::vector<chi_mesh::Vector3>& grad_shape_values,
    const std::vector<double>& phi,
    size_t num_nodes)
{
  chi_mesh::Vector3 gradphi_h;
  for (size_t j=0; j < num_nodes; ++j)
    gradphi_h += phi[j] * grad_shape_values[j];

  return gradphi_h;
}

//###################################################################
/**Obtains a field function interpolant of the angular flux-gradient.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  GetGradPsiH(
    const std::vector<chi_mesh::Vector3>& grad_shape_values,
    const std::vector<VecDbl>& phi,
    const chi_mesh::Vector3 omega,
    size_t num_nodes,
    const std::vector<EllEmIndices>& m_to_ell_em_map)
{
  const size_t num_moms = phi.size();

  const auto phi_theta   = OmegaToPhiThetaSafe(omega);
  const double& varphi   = phi_theta.first;
  const double& theta    = phi_theta.second;

  chi_mesh::Vector3 grad_psi_h(0.0, 0.0, 0.0);
  for (size_t m=0; m<num_moms; ++m)
  {
    const EllEmIndices& ell_em = m_to_ell_em_map[m];
    const unsigned int ell = ell_em.ell;
    const int           em = ell_em.m;

    const double Y_ell_em = chi_math::Ylm(ell, em, varphi, theta);
    const double factor = ((2.0*ell + 1.0)/4.0/M_PI);

    for (size_t j=0; j < num_nodes; ++j)
      grad_psi_h += factor * Y_ell_em * phi[m][j] * grad_shape_values[j];
  }

  return grad_psi_h;
}


//###################################################################
/***/
void mcpartra::ResidualSourceA::ExportCellResidualMoments()
{
  chi_log.Log() << "Exporting Cell Residual Moments";

  auto& rng              = ref_solver.rng0;
  auto& pwl              = ref_solver.pwl;
  auto& fv               = fv_sdm;
  auto& cell_geom_info   = *cell_geometry_info;

  const double FOUR_PI   = 4.0*M_PI;
  typedef chi_mesh::Vector3 Vec3;

  //======================================== Define quadrature
  chi_math::ProductQuadrature quadrature;
  quadrature.InitializeWithGLC(4,4);

  const auto& w_n = quadrature.weights;
  const auto& omega_n = quadrature.omegas;
  const auto& abscis_n = quadrature.abscissae;
  const size_t num_q_points = w_n.size();

  quadrature.MakeHarmonicIndices(1,3);

  const auto& m_to_ell_em_map = quadrature.GetMomentToHarmonicsIndexMap();
  const size_t num_moments = m_to_ell_em_map.size();

  //======================================== Define uk_man
  chi_math::UnknownManager uk_man_r;
  for (size_t g=0; g < num_groups; ++g)
    uk_man_r.AddUnknown(chi_math::UnknownType::VECTOR_N, num_moments);

  //======================================== Init vector
  const size_t num_dofs = fv->GetNumLocalDOFs(uk_man_r);
  std::vector<double> r_moments_local(num_dofs, 0.0);

  //======================================== Loop over cells
  //Predefine N_i and grad_N_i vectors to prevent
  //unnecessary memory reallocation.
  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  chi_log.Log() << "Checkpoint A";

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

      for (size_t n=0; n<num_q_points; ++n)
      {
        auto& omega = omega_n.at(n);
        auto& varphi = abscis_n.at(n).phi;
        auto& theta  = abscis_n.at(n).theta;
        double r_avg_omega_n = 0.0;

        size_t num_points = ref_solver.options.resid_src_integration_N_y;
        for (size_t i=0; i < num_points; ++i)
        {
          auto x_i   = GetRandomPositionInCell(rng, cell_geom_info[k]);

          cell_pwl_view->ShapeValues(x_i, shape_values);
          cell_pwl_view->GradShapeValues(x_i,grad_shape_values);

          double phi = GetPhiH(shape_values, nodal_phi, num_nodes);
          Vec3   grad_phi = GetGradPhiH(grad_shape_values, nodal_phi, num_nodes);

          r_avg_omega_n += (1.0 / FOUR_PI) * (Q - siga * phi - omega.Dot(grad_phi));
        }//for i
        r_avg_omega_n /= num_points;

        for (size_t m=0; m<num_moments; ++m)
        {
          auto& ell_em = m_to_ell_em_map[m];
          unsigned int ell = ell_em.ell;
          int           em = ell_em.m;

          uint64_t dof_map = fv->MapDOFLocal(cell,0,uk_man_r,g,m);

          r_moments_local.at(dof_map) += w_n[n] *
                                         std::fabs(r_avg_omega_n) *
                                         chi_math::Ylm(ell, em, varphi, theta);
        }
      }//for n
    }//for g

    auto progress = chi_misc_utils::
      PrintIterationProgress(cell.local_id, grid->local_cells.size());

    if (not progress.empty())
      chi_log.Log() << "    " << progress << "% Complete";
  }//for cell

  chi_log.Log() << "Checkpoint B";

  //============================================= Export interior source
  //                                              as FieldFunction
  auto fv_sd = std::dynamic_pointer_cast<SpatialDiscretization>(fv);
  auto R_ff = std::make_shared<chi_physics::FieldFunction>(
    "Zrmoms",                                     //Text name
    fv_sd,                                        //Spatial Discretization
    &r_moments_local,                             //Data
    uk_man_r,                                     //Nodal variable structure
    0, 0);                                        //Reference variable and component

  R_ff->ExportToVTKFV("Y_Rmoms","Zrmoms",true);

  chi_log.Log() << "Done exporting Cell Residual Moments";
}


