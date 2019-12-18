#include "solver_montecarlon.h"
#include "../Source/ResidualSource/mc_rmc2_source.h"

#include <ChiMesh/Raytrace/raytracing.h>

#include <chi_log.h>
extern ChiLog chi_log;

#include <ChiPhysics/chi_physics.h>
extern ChiPhysics chi_physics_handler;

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyRMC(
  chi_montecarlon::Particle &prtcl,
  chi_mesh::Vector pf,
  chi_mesh::RayDestinationInfo& ray_dest_info)
{
  //======================================== Get cell information
  auto cell = grid->cells[prtcl.cur_cell_ind];
  int cell_local_ind = cell->cell_local_id;
  auto cell_pwl_view = pwl_discretization->MapFeView(prtcl.cur_cell_ind);
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];
  auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double siga = xs->sigma_ag[prtcl.egrp];
  double sigt = xs->sigma_tg[prtcl.egrp];

  double tracklength = (pf - prtcl.pos).Norm();

  double avg_weight = 0.0;

  //--------------------------------- Get Src and mappings
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


  //--------------------------------- Lambda getting phi
  auto phi_s = [](std::vector<double>& N_in,
                  int dofs, int rmap,
                  chi_montecarlon::ResidualSource2* rsrc,
                  int egrp)
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
  };

  //--------------------------------- Lambda getting gradphi
  auto gradphi_s = [](std::vector<chi_mesh::Vector>& Grad_in,
                      int dofs, int rmap,
                      chi_montecarlon::ResidualSource2* rsrc,
                      int egrp)
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
  };

  //--------------------------------- Lambda computing q
  auto q_s = [](double siga,double phi, double Q,
                chi_mesh::Vector& omega,chi_mesh::Vector& nabla_phi)
  {
    double retval = Q - siga*phi - omega.Dot(nabla_phi);

    return retval;
  };

  //================================= Loop over segments
  chi_mesh::Vector p_i = prtcl.pos;
  chi_mesh::Vector p_f = prtcl.pos;

  //========================== Compute q_i
  cell_pwl_view->ShapeValues(p_i,N);
  cell_pwl_view->GradShapeValues(p_i,Grad);
  double phi_i = phi_s(N,cell_pwl_view->dofs,rmap,src,prtcl.egrp);
  auto gradphi_i = gradphi_s(Grad,cell_pwl_view->dofs,rmap,src,prtcl.egrp);

  double q_i = q_s(siga,phi_i,0.0,prtcl.dir,gradphi_i);
  double q_f = q_i;

  for (auto seg_len : segment_lengths) //segment_length
  {
//    int n_s = std::max(std::ceil(5*seg_len*sigt),4.0);
    int n_s = 2;
    double s_L = seg_len/n_s;

    for (int subdiv=0; subdiv<n_s; subdiv++)
    {
      p_i = p_f; //Only p_f gets changed below
      q_i = q_f;

      //========================== Compute q_f
      p_f = p_i + prtcl.dir*s_L;

      cell_pwl_view->ShapeValues(p_f,N);
      cell_pwl_view->GradShapeValues(p_f,Grad);
      double phi_f = phi_s(N,cell_pwl_view->dofs,rmap,src,prtcl.egrp);
      auto gradphi_f = gradphi_s(Grad,cell_pwl_view->dofs,rmap,src,prtcl.egrp);

      q_f = q_s(siga,phi_f,0.0,prtcl.dir,gradphi_f);

      //========================== Computing a and b
      double a = q_i;
      double b = (q_f-q_i)/s_L;

      //========================== Computing Iq
      double Iq  = (a/sigt)*(exp(sigt*s_L)-1.0);
      Iq += (b/sigt/sigt)*(sigt*s_L - 1.0)*exp(sigt*s_L);
      Iq -= (b/sigt/sigt)*(0.0      - 1.0)*1.0;

      //========================== Computing w_f
      double w_f = exp(-sigt*s_L)*prtcl.w + exp(-sigt*s_L)*Iq;

      double w_avg = 0.5*(w_f + prtcl.w);
      prtcl.w = w_f;
      avg_weight += s_L*w_avg;

      //========================== Contribute to tally

      //Get shape values at center of contrib
      auto p_c = (p_f+p_i)*0.5;
      cell_pwl_view->ShapeValues(p_c,N);

      for (int dof=0; dof<cell_pwl_view->dofs; dof++)
      {
        int ir = map + dof*num_grps*num_moms + num_grps*0 + prtcl.egrp;
        double pwl_tally_contrib = s_L*w_avg*N[dof];

        phi_pwl_tally[ir]     += pwl_tally_contrib;
        phi_pwl_tally_sqr[ir] += pwl_tally_contrib*pwl_tally_contrib;
      }//for dof

    }//for subdivision

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
  cell_pwl_view->ShapeValues(pf,N);
  double phi_neg = phi_s(N,cell_pwl_view->dofs,rmap,src,prtcl.egrp);

  double phi_pos = 0.0;
  int neighbor = ray_dest_info.destination_face_neighbor;
  if (neighbor>0)
  {
    auto adj_cell = grid->cells[neighbor];
    int adj_cell_local_ind = adj_cell->cell_local_id;
    auto adj_cell_pwl_view = pwl_discretization->MapFeView(neighbor);

    adj_cell_pwl_view->ShapeValues(pf,N);
    int armap =
      (*src->resid_ff->local_cell_dof_array_address)[adj_cell_local_ind];
    phi_pos = phi_s(N,adj_cell_pwl_view->dofs,armap,src,prtcl.egrp);

    prtcl.w += (phi_neg - phi_pos);
  }


}