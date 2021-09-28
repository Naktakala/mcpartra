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
  PopulateMaterialData(int mat_id, int group_g, MaterialData &mat_data)
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
                   int egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  double phi = 0.0;
  for (int dof=0; dof<dofs; dof++)
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
                       int egrp)
{
  auto& cell = grid->local_cells[cell_local_id];

  auto& sdm = resid_ff->spatial_discretization;

  if (sdm->type != chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
    throw std::invalid_argument(std::string(__PRETTY_FUNCTION__) +
                                " Invalid spatial discretization.");

  auto pwl_sdm = std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(sdm);

  auto& uk_man = resid_ff->unknown_manager;

  chi_mesh::Vector3 gradphi;
  for (int dof=0; dof<dofs; dof++)
  {
    int64_t ir = pwl_sdm->MapDOFLocal(cell,dof,uk_man,0,egrp);

    gradphi = gradphi + (Grad_in[dof]*(*resid_ff->field_vector_local)[ir]);
  }//for dof

  return gradphi;
}