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
double mcpartra::ResidualSourceA::
  GetResidualFFPhi(std::vector<double> &N_in,
                   size_t dofs,
                   uint64_t cell_local_id,
                   size_t egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  double phi = 0.0;
  for (size_t dof=0; dof<dofs; dof++)
  {
    int64_t ir = pwl_sdm->MapDOFLocal(cell,dof,uk_man,0,egrp);

    phi += (*resid_ff->field_vector_local)[ir]*N_in[dof];
  }//for dof

  return phi;
}

//###################################################################
/**Obtains a field function interpolant of the flux-gradient.*/
chi_mesh::Vector3 mcpartra::ResidualSourceA::
  GetResidualFFGradPhi(std::vector<chi_mesh::Vector3>& Grad_in,
                       size_t dofs,
                       uint64_t cell_local_id,
                       size_t egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  chi_mesh::Vector3 gradphi;
  for (size_t dof=0; dof<dofs; dof++)
  {
    int64_t ir = pwl_sdm->MapDOFLocal(cell,dof,uk_man,0,egrp);

    gradphi = gradphi + (Grad_in[dof]*(*resid_ff->field_vector_local)[ir]);
  }//for dof

  return gradphi;
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

  auto Isotropic_PDF = [TWO_PI](double mu) {return (1.0/TWO_PI);};
  auto PDF = [TWO_PI,a,b](double mu) {return (1.0/TWO_PI) * exp( a + b*mu );};

  //======================================== Find domain size
  double max_psi = 0.0;
  max_psi = std::max(max_psi, PDF(-1.0));
  max_psi = std::max(max_psi, PDF( 1.0));

  //======================================== Rejection sample pdf for mu
  bool rejected = true;
  double sampled_mu = 1.0;
  for (int i=0; i<10000; ++i)
  {
    double mu       = rng.Rand()*2.0 - 1.0;
    double random_y = rng.Rand()*max_psi;

    if (random_y < PDF(mu))
    {
      rejected = false;
      sampled_mu = mu;
    }

    if (not rejected) break;
  }
  if (rejected)
    std::cout << "Shit happened!\n";

  double weight_correction = Isotropic_PDF(sampled_mu) / PDF(sampled_mu);

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

  //Sample direction
  double theta    = acos(sampled_mu);
  double varphi   = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  return std::make_pair(R*ref_dir, weight_correction);
}
