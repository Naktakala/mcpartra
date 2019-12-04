#include "solver_montecarlon.h"

#include <ChiMesh/Cell/cell.h>

#include <FiniteVolume/CellViews/fv_slab.h>
#include <FiniteVolume/CellViews/fv_polygon.h>
#include <FiniteVolume/CellViews/fv_polyhedron.h>

#include <PiecewiseLinear/CellViews/pwl_slab.h>

#include <ChiMesh/Raytrace/raytracing.h>

#include "../Source/ResidualSource/mc_rmc2_source.h"

typedef unsigned long long TULL;

#include <chi_log.h>
#include <chi_mpi.h>

extern ChiLog chi_log;
extern ChiMPI chi_mpi;

#include <ChiPhysics/chi_physics.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

extern ChiPhysics chi_physics_handler;

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTally(
  chi_montecarlon::Particle &prtcl,
  chi_mesh::Vector pf)
{
  auto cell = grid->cells[prtcl.cur_cell_ind];
  int cell_local_ind = cell->cell_local_id;

  int ir = cell_local_ind*num_grps + prtcl.egrp;

  double tracklength = (pf - prtcl.pos).Norm();

  double tally_contrib = tracklength*prtcl.w;

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

  if (make_pwld)
  {
    segment_lengths.clear();
    segment_lengths.push_back(tracklength);
    chi_mesh::PopulateRaySegmentLengths(grid, cell,
                                        segment_lengths,
                                        prtcl.pos, pf,prtcl.dir);



    auto cell_pwl_view = pwl_discretization->MapFeView(prtcl.cur_cell_ind);
    int map            = local_cell_pwl_dof_array_address[cell_local_ind];

    double last_segment_length = 0.0;

    for (auto segment_length : segment_lengths)
    {
      double d = last_segment_length + 0.5*segment_length;
      last_segment_length += segment_length;
      auto p = prtcl.pos + prtcl.dir*d;

      cell_pwl_view->ShapeValues(p,N);

      for (int dof=0; dof<cell_pwl_view->dofs; dof++)
      {
        ir = map + dof*num_grps*num_moms + num_grps*0 + prtcl.egrp;
        double pwl_tally_contrib = segment_length*prtcl.w*N[dof];

        phi_pwl_tally[ir]     += pwl_tally_contrib;
        phi_pwl_tally_sqr[ir] += pwl_tally_contrib*pwl_tally_contrib;
      }//for dof
    }//for segment_length
  }//if make pwld

}

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyRMC(
  chi_montecarlon::Particle &prtcl,
  chi_mesh::Vector pf)
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

  //======================================== Contribute PWLD
  if (make_pwld)
  {
    //--------------------------------- Develop segments
    segment_lengths.clear();
    segment_lengths.push_back(tracklength);
    chi_mesh::PopulateRaySegmentLengths(grid, cell,
                                        segment_lengths,
                                        prtcl.pos, pf,prtcl.dir);

    //--------------------------------- Get Src and mappings
    auto src = (chi_montecarlon::ResidualSource2*)sources.back();
    int map  = local_cell_pwl_dof_array_address[cell_local_ind];
    int rmap = (*src->resid_ff->local_cell_dof_array_address)[cell_local_ind];

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
      const double four_pi = 4.0*acos(-1.0);

      double retval = Q - siga*phi - omega.Dot(nabla_phi);

      return retval/four_pi;
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
      int n_s = std::max(std::ceil(5*seg_len*sigt),4.0);
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

//        if (prtcl.dir.z <0.0)
//        {
//          chi_log.Log(LOG_0)
//            << " pos_i=" << p_i.PrintS()
//            << " pos_f=" << p_f.PrintS()
//            << " phi_i=" << phi_i
//            << " phi_f=" << phi_f
//            << " q_i=" << q_i
//            << " q_f=" << q_f
//            << " w=" << prtcl.w;
//
//          usleep(1000000);
//        }


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

//          if (prtcl.dir.z <0.0)
//            chi_log.Log(LOG_0) << w_avg << " " << N[dof];
        }//for dof

//        if (prtcl.dir.z <0.0)
//          usleep(100000);
      }//for subdivision

    }//for segment_length
  }//if make pwld
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

}

//###################################################################
/**Merges tallies from multiple locations.*/
void chi_montecarlon::Solver::RendesvouzTallies()
{
  TULL temp_nps_global = 0;
  MPI_Allreduce(&nps,&temp_nps_global,1,MPI_UNSIGNED_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);

  nps_global += temp_nps_global;

  //============================================= If mesh global
  if (mesh_is_global)
  {
    //============================ Merge phi_local
    int num_values = phi_tally.size();
    std::vector<double> temp_phi_global(num_values,0.0);

    MPI_Allreduce(phi_tally.data(),
                  temp_phi_global.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    for (int i=0; i<num_values; i++)
      phi_global[i] += temp_phi_global[i];

    //============================ Merge the square of phi_local
    phi_global_tally_sqr.assign(num_values,0.0);

    MPI_Allreduce(phi_tally_sqr.data(),
                  phi_global_tally_sqr.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    //============================ Reset tallies
    phi_tally.assign(phi_tally.size(),0.0);
    phi_tally_sqr.assign(phi_tally_sqr.size(),0.0);
  }
  //============================================= Mesh partitioned
  else
  {
    int num_values = phi_tally.size();

    phi_global_tally_sqr.assign(num_values,0.0);

    for (int i=0; i<num_values; i++)
    {
      phi_global[i] += phi_tally[i];
      phi_global_tally_sqr[i] = phi_tally_sqr[i];
    }

    //============================ Reset tallies
    phi_tally.assign(phi_tally.size(),0.0);
    phi_tally_sqr.assign(phi_tally_sqr.size(),0.0);
  }


}

//###################################################################
/**Merges tallies from multiple locations.*/
void chi_montecarlon::Solver::RendesvouzPWLTallies()
{
  //============================================= If mesh global
  if (mesh_is_global)
  {
    //============================ Merge phi_local
    int num_values = phi_pwl_tally.size();
    std::vector<double> temp_phi_pwl_global(num_values,0.0);

    MPI_Allreduce(phi_pwl_tally.data(),
                  temp_phi_pwl_global.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    for (int i=0; i<num_values; i++)
      phi_pwl_global[i] += temp_phi_pwl_global[i];

    //============================ Merge square if phi_local
    phi_pwl_global_tally_sqr.assign(num_values,0.0);

    MPI_Allreduce(phi_pwl_tally_sqr.data(),
                  phi_pwl_global_tally_sqr.data(),
                  num_values,MPI_DOUBLE,
                  MPI_SUM,MPI_COMM_WORLD);

    //============================ Reset tallies
    phi_pwl_tally.assign(phi_pwl_tally.size(),0.0);
    phi_pwl_tally_sqr.assign(phi_pwl_tally_sqr.size(),0.0);
  }
  //============================================= Mesh partitioned
  else
  {
    int num_values = phi_pwl_tally.size();

    phi_pwl_global_tally_sqr.assign(num_values,0.0);

    for (int i=0; i<num_values; i++)
    {
      phi_pwl_global[i] += phi_pwl_tally[i];
      phi_pwl_global_tally_sqr[i] = phi_pwl_tally_sqr[i];
    }

    //============================ Reset tallies
    phi_pwl_tally.assign(phi_pwl_tally.size(),0.0);
    phi_pwl_tally_sqr.assign(phi_pwl_tally_sqr.size(),0.0);
  }
}


//###################################################################
/**Computes the relative std dev for all the tallies.*/
void chi_montecarlon::Solver::ComputeRelativeStdDev()
{
  max_relative_error = 0.0;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=0; g<num_grps; g++)
    {
      int ir = lc*num_grps + g;

      TULL n = batch_sizes[current_batch];

      double x_avg  = phi_global[ir]/nps_global;
      double x2_avg = (phi_global_tally_sqr[ir]/n);

      double stddev = sqrt((x2_avg - x_avg*x_avg)/nps_global);

      if (!std::isinf(stddev/x_avg))
      {
        phi_local_relsigma[ir] = stddev/x_avg;

        if (g==0 && lc == 0)
        {
          max_relative_error = phi_local_relsigma[ir];
          max_relative_error2 = x2_avg;
          max_relative_error3 = x_avg*x_avg;
        }
      }
      else
      {
        phi_local_relsigma[ir] = 0.0;
      }


      if (std::isinf(phi_local_relsigma[ir]))
      {
        printf("Infinite lc=%d g=%d\n", lc,g);
        chi_log.Log(LOG_ALL) << "stddev=" << stddev << "\n";
        chi_log.Log(LOG_ALL) << "xavg=" << x_avg << "\n";
        chi_log.Log(LOG_ALL) << "x2avg=" << x2_avg << "\n";
        chi_log.Log(LOG_ALL) << "phirel=" << phi_local_relsigma[ir] << "\n";
        chi_log.Log(LOG_ALL) << "batch=" << n << "\n";
        chi_log.Log(LOG_ALL) << "nps=" << nps_global << "\n";
      }
    }//for g
  }//for local cell
}


//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::NormalizeTallies()
{
  if (nps_global == 0) nps_global = 1;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell = grid->cells[cell_glob_index];

    double V = 1.0;

    if (cell->Type() == chi_mesh::CellType::SLAB)
    {
      auto cell_fv_view =
        (SlabFVView*)fv_discretization->MapFeView(cell->cell_global_id);
      V = cell_fv_view->volume;
    }
    else if (cell->Type() == chi_mesh::CellType::POLYGON)
    {
      auto cell_fv_view =
        (PolygonFVView*)fv_discretization->MapFeView(cell->cell_global_id);
      V = cell_fv_view->volume;
    }
    else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
    {
      auto cell_fv_view =
        (PolyhedronFVView*)fv_discretization->MapFeView(cell->cell_global_id);
      V = cell_fv_view->volume;
    }
    else
    {
      chi_log.Log(LOG_ALLERROR)
        << "Unsupported cell type encountered in call to "
        << "chi_montecarlon::Solver::ComputeTallySqr.";
      exit(EXIT_FAILURE);
    }


    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=hi; g<=lo; g++)
    {
      int ir = lc*num_grps + g;

      phi_global[ir] *= tally_multipl_factor/nps_global/V;
      phi_global[ir] += phi_global_initial_value[ir];
      if (not phi_uncollided_rmc.empty())
        phi_global[ir] += phi_uncollided_rmc[ir];
    }//for g
  }//for local cell
}

//###################################################################
/**Compute tally square contributions.*/
void chi_montecarlon::Solver::NormalizePWLTallies()
{
  if (nps_global == 0) nps_global = 1;

  int num_cells = grid->local_cell_glob_indices.size();
  for (int lc=0; lc<num_cells; lc++)
  {
    int cell_glob_index = grid->local_cell_glob_indices[lc];
    auto cell_pwl_view   = pwl_discretization->MapFeView(cell_glob_index);
    int map             = local_cell_pwl_dof_array_address[lc];

    int hi = 0;
    int lo = num_grps-1;
    if (group_hi_bound >= 0)
      hi = group_hi_bound;
    if (group_lo_bound >=0 && group_lo_bound >= group_hi_bound)
      lo = group_lo_bound;

    for (int g=hi; g<=lo; g++)
    {
      for (int dof=0; dof<cell_pwl_view->dofs; dof++)
      {
        int ir = map + dof*num_grps*num_moms + num_grps*0 + g;
        double V = cell_pwl_view->IntV_shapeI[dof];

        phi_pwl_global[ir] *= tally_multipl_factor/nps_global/V;

      }//for dof
    }//for g
  }//for local cell
}