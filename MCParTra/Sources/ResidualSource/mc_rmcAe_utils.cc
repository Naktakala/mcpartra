#include "mc_rmcA_source.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <ChiPhysics/chi_physics.h>
#include "chi_log.h"
#include "SourceDrivenSolver/sdsolver.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Populates a material information data structure from a mat-id.*/
void mcpartra::ResidualSourceA::
  PopulateMaterialData(int mat_id, size_t group_g, MaterialData &mat_data)
{
  auto xs = ref_solver.matid_xs_map2[mat_id];
  double siga = xs->sigma_a[group_g];

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
/**Samples and exponential importance representation.*/
std::pair<chi_mesh::Vector3,double> mcpartra::ResidualSourceA::
  SampleSpecialRandomDirection(chi_math::RandomNumberGenerator &rng,
                               size_t group,
                               uint64_t cell_local_id)
{
  size_t dof_map = cell_local_id*num_groups + group;

  const auto& omega_J  = ref_solver.local_cell_importance_directions[dof_map];
  const auto& a_b_pair = ref_solver.local_cell_importance_exp_coeffs[dof_map];
  const double a = a_b_pair.first;
  const double b = a_b_pair.second;

  if (std::fabs(b) < 1.0e-8)
    return std::make_pair(SampleRandomDirection(rng),1.0);

  //======================================== Define utilities
  const double TWO_PI = 2.0*M_PI;

  auto Isotropic_PDF = [](double mu) {return 0.5;};
  auto PDF = [TWO_PI,a,b](double mu) {return TWO_PI * exp( a + b*mu );};

  //======================================== Rejection sample pdf for mu
  // Find domain size
  double max_psi = 0.0;
  max_psi = std::max(max_psi, PDF(-1.0));
  max_psi = std::max(max_psi, PDF( 1.0));

  bool rejected = true;
  double mu_prime = 1.0;
  for (int i=0; i<10000; ++i)
  {
    mu_prime = rng.Rand() * 2.0 - 1.0;
    const double random_PDF = rng.Rand() * max_psi;

    if (random_PDF < PDF(mu_prime)) rejected = false;
    if (not rejected) break;
  }

//  const double C_0 = (TWO_PI/b) * exp(a - b);
//  const double theta_dvi_C0 = rng.Rand()/C_0;
//  double mu = (1.0/b) * log( exp(-b) * (theta_dvi_C0 + 1) );

  double weight_correction = Isotropic_PDF(mu_prime) / PDF(mu_prime);

  //======================================== Compute omega in ref-coordinates
  //Sample direction
  double theta  = acos(mu_prime);
  double varphi = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 omega_prime;
  omega_prime.x = sin(theta) * cos(varphi);
  omega_prime.y = sin(theta) * sin(varphi);
  omega_prime.z = cos(theta);

  //======================================== Perform rotation
  //Build rotation matrix
  chi_mesh::Matrix3x3 R;

  const chi_mesh::Vector3 khat(0.0,0.0,1.0);

  if      (omega_J.Dot(khat) >  0.9999999)
    R.SetDiagonalVec(1.0, 1.0, 1.0);
  else if (omega_J.Dot(khat) < -0.9999999)
    R.SetDiagonalVec(1.0, 1.0,-1.0);
  else
  {
    chi_mesh::Vector3 binorm = khat.Cross(omega_J);
    binorm = binorm/binorm.Norm();

    chi_mesh::Vector3 tangent = binorm.Cross(omega_J);
    tangent = tangent/tangent.Norm();

    R.SetColJVec(0, tangent);
    R.SetColJVec(1, binorm);
    R.SetColJVec(2, omega_J);
  }

  chi_mesh::Vector3 omega = R * omega_prime;

  return std::make_pair(omega, weight_correction);
}
