#include "solver_montecarlon.h"

#include <ChiMesh/Cell/cell.h>
#include <ChiMesh/Raytrace/raytracing.h>

#include <chi_log.h>
extern ChiLog chi_log;

#include <ChiPhysics/chi_physics.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

extern ChiPhysics chi_physics_handler;

typedef unsigned long long TULL;

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

