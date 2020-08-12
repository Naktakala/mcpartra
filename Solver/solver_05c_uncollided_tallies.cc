#include "solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyUNC(
  chi_montecarlon::Particle &prtcl,
  chi_mesh::Vector3 pf,
  double sig_t)
{
  auto cell = &grid->local_cells[prtcl.cur_cell_local_id];
  int cell_local_ind = cell->local_id;

  int ir = cell_local_ind*num_grps + prtcl.egrp;

  double tracklength = (pf - prtcl.pos).Norm();

  double avg_w = (sig_t<1.0e-16)? prtcl.w :
                 prtcl.w*(1.0 - exp(-sig_t*tracklength))/sig_t/tracklength;

  double w_exit = prtcl.w*exp(-sig_t*tracklength);

  double tally_contrib = tracklength*avg_w;

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



    auto cell_pwl_view = pwl_discretization->MapFeViewL(cell_local_ind);
    int map            = local_cell_pwl_dof_array_address[cell_local_ind];

    cell_pwl_view->ShapeValues(prtcl.pos, N_i);

    double last_segment_length = 0.0;
    for (auto segment_length : segment_lengths)
    {
      double d = last_segment_length + segment_length;
      last_segment_length += segment_length;
      auto p = prtcl.pos + prtcl.dir*d;

      cell_pwl_view->ShapeValues(p, N_f);

      int moment = 0;
      for (int dof=0; dof<cell_pwl_view->dofs; dof++)
      {
        ir = map + dof*num_grps*num_moms + num_grps*moment + prtcl.egrp;

        double ell = segment_length;

        double w_avg  = (N_i[dof]/sig_t)*(1.0-exp(-sig_t*ell));
               w_avg += ((N_f[dof]-N_i[dof])/(sig_t*sig_t*ell))*
                 (1.0 - (1+sig_t*ell)*exp(-sig_t*ell));
               w_avg *= prtcl.w/ell;

        double pwl_tally_contrib = segment_length * w_avg;

        phi_pwl_tally[ir]     += pwl_tally_contrib;
        phi_pwl_tally_sqr[ir] += pwl_tally_contrib*pwl_tally_contrib;
      }//for dof

      //reset for new segment
      N_i = N_f;
      prtcl.w *= exp(-sig_t*segment_length);
    }//for segment_length
  }//if make pwld
  else
    prtcl.w = w_exit;
}