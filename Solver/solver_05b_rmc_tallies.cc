#include "solver_montecarlon.h"
#include "../Source/ResidualSource/mc_rmc2_source.h"

#include <ChiMesh/Raytrace/raytracing.h>
#include "ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h"

#include <chi_log.h>
extern ChiLog chi_log;

#include <ChiPhysics/chi_physics.h>
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Obtains a field function interpolant of the flux.*/
double chi_montecarlon::Solver::
  GetResidualFFPhi(std::vector<double> &N_in, int dofs, int rmap,
                   chi_montecarlon::ResidualSource2 *rsrc, int egrp)
{
  double phi = 0.0;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = rmap +
             dof*rsrc->resid_ff->num_grps*rsrc->resid_ff->num_moms +
             rsrc->resid_ff->num_grps*0 +
             egrp;
    phi += (*rsrc->resid_ff->field_vector_local)[ir]*N_in[dof];
  }//for dof

  return phi;
}

//###################################################################
/**Obtains a field function interpolant of the flux-gradient.*/
chi_mesh::Vector chi_montecarlon::Solver::
GetResidualFFGradPhi(std::vector<chi_mesh::Vector>& Grad_in, int dofs, int rmap,
                 chi_montecarlon::ResidualSource2 *rsrc, int egrp)
{
  chi_mesh::Vector gradphi;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = rmap +
             dof*rsrc->resid_ff->num_grps*rsrc->resid_ff->num_moms +
             rsrc->resid_ff->num_grps*0 +
             egrp;
    gradphi = gradphi + (Grad_in[dof]*(*rsrc->resid_ff->field_vector_local)[ir]);
  }//for dof

  return gradphi;
}


//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyRMC(
  chi_montecarlon::Particle &prtcl,
  chi_mesh::Vector pf,
  chi_mesh::RayDestinationInfo& ray_dest_info)
{
  //======================================== Get cell information
  auto cell           = grid->cells[prtcl.cur_cell_ind];
  int  cell_local_ind = cell->cell_local_id;
  auto cell_pwl_view  = pwl_discretization->MapFeView(prtcl.cur_cell_ind);

  //======================================== Get material properties
  int  mat_id         = cell->material_id;
  int  xs_prop_id     = matid_xs_map[mat_id];
  int  src_prop_id    = matid_q_map[mat_id];
  auto material = chi_physics_handler.material_stack[mat_id];
  auto xs = (chi_physics::TransportCrossSections*)material->properties[xs_prop_id];

  double siga = xs->sigma_ag[prtcl.egrp];
  double sigt = xs->sigma_tg[prtcl.egrp];
  double Q    = 0.0;
  if (src_prop_id >= 0)
  {
    auto prop = material->properties[src_prop_id];
    auto q_prop = (chi_physics::IsotropicMultiGrpSource*)prop;
    Q = q_prop->source_value_g[prtcl.egrp];
  }

  //======================================== Compute tracklength and init weight
  double tracklength = (pf - prtcl.pos).Norm();
  double avg_weight = 0.0;

  //======================================== Get residual-source and mappings
  auto src = (chi_montecarlon::ResidualSource2*)sources.back();
  int map  = local_cell_pwl_dof_array_address[cell_local_ind];
  int rmap = (*src->resid_ff->local_cell_dof_array_address)[cell_local_ind];

  //======================================== Contribute PWLD
  //--------------------------------- Develop segments
  segment_lengths.clear();
  segment_lengths.push_back(tracklength);
  chi_mesh::PopulateRaySegmentLengths(grid, cell,
                                      segment_lengths,
                                      prtcl.pos, pf,prtcl.dir);

  //--------------------------------- Lambda computing residual source
  auto Resdiual_Q = [](double siga, double phi, double Q,
                       chi_mesh::Vector& omega, chi_mesh::Vector& nabla_phi)
  {
    double retval = Q - siga*phi - omega.Dot(nabla_phi);

    return retval;
  };

  //--------------------------------- Pre-fetch shape-function values
  chi_mesh::Vector p_i = prtcl.pos;
  chi_mesh::Vector p_f = prtcl.pos;

  cell_pwl_view->ShapeValues(p_i, N_f);
  cell_pwl_view->GradShapeValues(p_i,Grad);

  //--------------------------------- Compute q_i
  double phi_i = GetResidualFFPhi(N_f, cell_pwl_view->dofs, rmap, src, prtcl.egrp);
  auto gradphi_i = GetResidualFFGradPhi(Grad,cell_pwl_view->dofs,rmap,src,prtcl.egrp);

  //Compute residual source at i
  double q_i = Resdiual_Q(siga, phi_i, Q, prtcl.dir, gradphi_i);
  double q_f = q_i;

  //--------------------------------- Loop over segments
  //Each segment starts at p_i and ends at p_f
  //The general notation is _i = begin, _f = end
  for (auto s_L : segment_lengths) //segment_length
  {
    p_i = p_f; //Only p_f gets changed below
    q_i = q_f;
    N_i = N_f;

    //========================== Compute q_f
    p_f = p_i + prtcl.dir*s_L;

    //Grab shape function values at segment end
    cell_pwl_view->ShapeValues(p_f, N_f);
    cell_pwl_view->GradShapeValues(p_f,Grad);

    //Grab phi and grad-phi from residual's
    //reference field function
    double phi_f = GetResidualFFPhi(N_f, cell_pwl_view->dofs,
                                    rmap, src, prtcl.egrp);
    auto gradphi_f = GetResidualFFGradPhi(Grad,cell_pwl_view->dofs,
                                          rmap,src,prtcl.egrp);

    //Compute residual source at f
    q_f = Resdiual_Q(siga, phi_f, Q, prtcl.dir, gradphi_f);

    //========================== Computing equation coefficients
    double w_i = prtcl.w;
    double a = q_i;
    double b = (q_f-q_i)/s_L;

    double c_0 = a/sigt - b/(sigt*sigt);
    double c_1 = b/sigt;
    double c_2 = w_i - c_0;

    //========================== Computing Iq
    double Iq  = (a/sigt)*(exp(sigt*s_L)-1.0);
    Iq += (b/sigt/sigt)*(sigt*s_L - 1.0)*exp(sigt*s_L);
    Iq -= (b/sigt/sigt)*(0.0      - 1.0)*1.0;

    //========================== Computing w_f
    double w_f = exp(-sigt*s_L)*w_i + exp(-sigt*s_L)*Iq;
    prtcl.w = w_f;

    //========================== Contribute to tally and compute average
    double w_avg = 0.0;
    for (int i=0; i < cell_pwl_view->dofs; ++i)
    {
      double c_3 = (N_f[i] - N_i[i]) / s_L;

      double c_4 = c_0 * N_i[i] * s_L;
      double c_5 = (c_0*c_3 + c_1*N_i[i])*s_L*s_L/2.0;
      double c_6 = c_1*c_3*s_L*s_L*s_L/3.0;
      double c_7 = -c_2*N_i[i]*(exp(-sigt*s_L)-1.0)/sigt;
      double c_8 = c_2*c_3/(sigt*sigt);
      double c_9 = 1.0 - (s_L*sigt+1.0)*exp(-sigt*s_L);

      double w_avg_i = (c_4+c_5+c_6+c_7+c_8*c_9)/s_L;

      int ir = map + i * num_grps * num_moms + num_grps * 0 + prtcl.egrp;
      double pwl_tally_contrib = s_L * w_avg_i;

      phi_pwl_tally[ir]     += pwl_tally_contrib;
      phi_pwl_tally_sqr[ir] += pwl_tally_contrib*pwl_tally_contrib;

      w_avg += w_avg_i;
    }//for dof

    avg_weight += s_L*w_avg;
  }//for segment_length
  avg_weight/=tracklength;

  //======================================== Contribute avg tally
  int ir = cell_local_ind*num_grps + prtcl.egrp;

  double tally_contrib = tracklength*avg_weight;

  phi_tally[ir]     += tally_contrib;
  phi_tally_sqr[ir] += tally_contrib*tally_contrib;

  if (std::isnan(tracklength))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Tracklength corruption."
      << " pos  " << prtcl.pos.PrintS()
      << " posf " << pf.PrintS();
    exit(EXIT_FAILURE);
  }

  //======================================== Adjust weight for face jump
  cell_pwl_view->ShapeValues(pf, N_f);
  double phi_neg = GetResidualFFPhi(N_f, cell_pwl_view->dofs, rmap, src, prtcl.egrp);

  double phi_pos = 0.0;
  int neighbor = ray_dest_info.destination_face_neighbor;
  if (neighbor>0)
  {
    auto adj_cell = grid->cells[neighbor];
    int adj_cell_local_ind = adj_cell->cell_local_id;
    auto adj_cell_pwl_view = pwl_discretization->MapFeView(neighbor);

    adj_cell_pwl_view->ShapeValues(pf, N_f);
    int armap =
      (*src->resid_ff->local_cell_dof_array_address)[adj_cell_local_ind];
    phi_pos = GetResidualFFPhi(N_f, adj_cell_pwl_view->dofs, armap, src, prtcl.egrp);

    prtcl.w += (phi_neg - phi_pos);
  }


}