#include "solver_montecarlon.h"
#include "../Source/ResidualSource/mc_rmcB_source.h"

#include <ChiMesh/Raytrace/raytracing.h>
#include "ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h"

#include <chi_log.h>
extern ChiLog& chi_log;

#include <ChiPhysics/chi_physics.h>
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Obtains a field function interpolant of the flux.*/
double chi_montecarlon::Solver::
  GetResidualFFPhi(std::vector<double> &N_in, int dofs, int rmap,
                   chi_montecarlon::ResidualSourceB *rsrc, int egrp)
{
  double phi = 0.0;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = rmap +
             dof * rsrc->resid_ff->num_components * rsrc->resid_ff->num_sets +
             rsrc->resid_ff->num_components * 0 +
             egrp;
    phi += (*rsrc->resid_ff->field_vector_local)[ir]*N_in[dof];
  }//for dof

  return phi;
}

//###################################################################
/**Obtains a field function interpolant of the flux-gradient.*/
chi_mesh::Vector3 chi_montecarlon::Solver::
GetResidualFFGradPhi(std::vector<chi_mesh::Vector3>& Grad_in, int dofs, int rmap,
                     chi_montecarlon::ResidualSourceB *rsrc, int egrp)
{
  chi_mesh::Vector3 gradphi;
  for (int dof=0; dof<dofs; dof++)
  {
    int ir = rmap +
             dof * rsrc->resid_ff->num_components * rsrc->resid_ff->num_sets +
             rsrc->resid_ff->num_components * 0 +
             egrp;
    gradphi = gradphi + (Grad_in[dof]*(*rsrc->resid_ff->field_vector_local)[ir]);
  }//for dof

  return gradphi;
}


//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyRMC(
  chi_montecarlon::Particle &prtcl,
  const chi_mesh::Vector3& pf,
  chi_mesh::RayDestinationInfo& ray_dest_info)
{
  //======================================== Get cell information
  auto cell           = &grid->local_cells[prtcl.cur_cell_local_id];
  int  cell_local_ind = cell->local_id;
  auto cell_pwl_view  = pwl->MapFeViewL(cell_local_ind);

  //======================================== Get material properties
  int  mat_id         = cell->material_id;
  int  xs_prop_id     = matid_xs_map[mat_id];
  int  src_prop_id    = matid_q_map[mat_id];
  auto material = chi_physics_handler.material_stack[mat_id];
  auto xs = (chi_physics::TransportCrossSections*)material->properties[xs_prop_id];

  double siga = xs->sigma_ag[prtcl.egrp];
  double sigt = xs->sigma_tg[prtcl.egrp];
  double sigs = sigt-siga;
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
  auto src = (chi_montecarlon::ResidualSourceB*)sources.back();
  int rmap = (*src->resid_ff->local_cell_dof_array_address)[cell_local_ind];

  //======================================== Contribute PWLD
  //--------------------------------- Develop segments
  segment_lengths.clear();
  segment_lengths.push_back(tracklength);
  chi_mesh::PopulateRaySegmentLengths(*grid, *cell,
                                      segment_lengths,
                                      prtcl.pos, pf,prtcl.dir);

  //--------------------------------- Lambda computing residual source
  auto Resdiual_Q = [](double siga, double phi, double Q,
                       chi_mesh::Vector3& omega, chi_mesh::Vector3& nabla_phi)
  {
    double retval = Q - siga*phi - omega.Dot(nabla_phi);

    return retval;
  };

  auto Sign = [](double a)
  {
    if (a<0.0) return -1.0;
    else return 1.0;
  };

  //--------------------------------- Pre-fetch shape-function values
  chi_mesh::Vector3 p_i = prtcl.pos;
  chi_mesh::Vector3 p_f = prtcl.pos;

  cell_pwl_view->ShapeValues(p_i, N_f);
  cell_pwl_view->GradShapeValues(p_i,Grad);

  //--------------------------------- Compute q_i
  double phi_i = GetResidualFFPhi(N_f, cell_pwl_view->dofs, rmap, src, prtcl.egrp);
  auto gradphi_i = GetResidualFFGradPhi(Grad,cell_pwl_view->dofs,rmap,src,prtcl.egrp);

  if (prtcl.pre_cell_global_id >= 0)
    prtcl.w -= phi_i;

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
    double Iq  = (a/sigt)*(exp(std::fmin(sigt*s_L,100.0))-1.0);
    Iq += (b/sigt/sigt)*(sigt*s_L - 1.0)*exp(std::fmin(sigt*s_L,100.0));
    Iq -= (b/sigt/sigt)*(0.0      - 1.0)*1.0;

    //========================== Computing w_f
    double w_f_before = prtcl.w;
    double w_f = exp(-sigt*s_L)*w_i + exp(-sigt*s_L)*Iq;
    prtcl.w = w_f;

    if (std::isnan(w_f) or std::isinf(w_f))
    {
      chi_log.Log(LOG_ALLERROR)
        << "Segfault in w, before = " << w_f_before << "\n"
        << "Iq=" << Iq << "\n"
        << "s_L=" << s_L;
      exit(EXIT_FAILURE);
    }

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

      int ir = pwl->MapDFEMDOFLocal(cell,i,&uk_man_fem,/*m*/0,prtcl.egrp);
      double pwl_tally_contrib = s_L * w_avg_i;

      for (int t : pwl_tallies)
      {
        if (not (prtcl.tally_mask & (1 << t))) continue;

        grid_tally_blocks[t].tally_local[ir]     += pwl_tally_contrib;
        grid_tally_blocks[t].tally_sqr_local[ir] += pwl_tally_contrib*pwl_tally_contrib;
      }

      w_avg += w_avg_i;
    }//for dof

    avg_weight += s_L*w_avg;
  }//for segment_length
  avg_weight/=tracklength;

  //======================================== Contribute avg tally
  int ir = fv->MapDOFLocal(cell,&uk_man_fv,/*m*/0,prtcl.egrp);

  double tally_contrib = tracklength*avg_weight;

  for (int t : fv_tallies)
  {
    if (not (prtcl.tally_mask & (1 << t))) continue;

    grid_tally_blocks[t].tally_local[ir]     += tally_contrib;
    grid_tally_blocks[t].tally_sqr_local[ir] += tally_contrib*tally_contrib;
  }

  if (prtcl.tally_mask & TallyMask::MAKE_DIRECT_PARTICLES)
  {
    double q_abs = std::fabs(tracklength*avg_weight*sigs);

    if (q_abs < 1.0)
    {
      if (rng0.Rand() < q_abs)
        particle_source_bank.push_back(
          MakeScatteredParticle(prtcl,tracklength,Sign(avg_weight)));
    }
    else if (q_abs >= 1.0)
    {
      int definite_amount = std::floor(q_abs);
      double fractional_amount = q_abs - definite_amount;

      for (int p=0; p< definite_amount; ++p)
      {
        particle_source_bank.push_back(
          MakeScatteredParticle(prtcl,tracklength,Sign(avg_weight)));
      }
      if (rng0.Rand() < fractional_amount)
        particle_source_bank.push_back(
          MakeScatteredParticle(prtcl,tracklength,Sign(avg_weight)));
    }
  }

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

  prtcl.w += phi_neg;
}

//###################################################################
/**Makes a scattered particle.*/
chi_montecarlon::Particle chi_montecarlon::Solver::
  MakeScatteredParticle(Particle &prtcl,
                        const double tracklength,
                        const double weight)
{
  Particle new_particle;

  //======================================== Sample position
  new_particle.pos = prtcl.pos + rng0.Rand()*tracklength*prtcl.dir;

  //======================================== Sample direction
  double costheta = 2.0*rng0.Rand()-1.0;
  double theta    = acos(costheta);
  double varphi   = rng0.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;


  //======================================== Weight, Energy, Liveliness
  new_particle.w = weight;
  new_particle.egrp = prtcl.egrp;

  new_particle.alive = true;

  if (std::fabs(weight) < 1.0e-12)
    new_particle.alive = false;

  //======================================== Cell ownership
  new_particle.cur_cell_global_id = prtcl.cur_cell_global_id;
  new_particle.cur_cell_local_id  = prtcl.cur_cell_local_id;

  //======================================== Methods
  new_particle.ray_trace_method = chi_montecarlon::Solver::RayTraceMethod::STANDARD;
  new_particle.tally_method = chi_montecarlon::Solver::TallyMethod::STANDARD;

  return new_particle;
}