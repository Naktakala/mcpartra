#include "solver_montecarlon.h"

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Makes a contribution to tallies*/
void chi_montecarlon::Solver::ContributeTallyUNC(
  chi_montecarlon::Particle &prtcl,
  const chi_mesh::Vector3& pf,
  double sig_t)
{
  auto cell = &grid->local_cells[prtcl.cur_cell_local_id];
  int cell_local_ind = cell->local_id;

  int ir = fv->MapDOFLocal(cell,&uk_man_fv,/*m*/0,prtcl.egrp);

  double tracklength = (pf - prtcl.pos).Norm();

  double avg_w = (sig_t<1.0e-16)? prtcl.w :
                 prtcl.w*(1.0 - exp(-sig_t*tracklength))/sig_t/tracklength;

  double w_exit = prtcl.w*exp(-sig_t*tracklength);

  double tally_contrib = tracklength*avg_w;

  //============================================= FV Tallies
  for (int t : fv_tallies)
  {
    if (prtcl.tally_mask & (1 << t))
    {
      grid_tally_blocks[t].tally_local[ir]     += tally_contrib;
      grid_tally_blocks[t].tally_sqr_local[ir] += tally_contrib*
                                                  tally_contrib;
    }//if tally applies
  }//for fv tallies

  if (std::isnan(tracklength))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Tracklength corruption."
      << " pos  " << prtcl.pos.PrintS()
      << " posf " << pf.PrintS();
    exit(EXIT_FAILURE);
  }

  //============================================= PWL Tallies
  for (int t : pwl_tallies)
  {
    if ( (prtcl.tally_mask & (1 << t)) && (make_pwld))
    {
      segment_lengths.clear();
      segment_lengths.push_back(tracklength);
      chi_mesh::PopulateRaySegmentLengths(*grid, *cell,
                                          segment_lengths,
                                          prtcl.pos, pf,prtcl.dir);

      auto cell_pwl_view = pwl->MapFeViewL(cell_local_ind);

      cell_pwl_view->ShapeValues(prtcl.pos, N_i);

      double last_segment_length = 0.0;
      for (auto segment_length : segment_lengths)
      {
        double d = last_segment_length + segment_length;
        last_segment_length += segment_length;
        auto p = prtcl.pos + prtcl.dir*d;

        cell_pwl_view->ShapeValues(p, N_f);

        for (int dof=0; dof<cell_pwl_view->dofs; dof++)
        {
          ir = pwl->MapDFEMDOFLocal(cell,dof,&uk_man_fem,/*m*/0,prtcl.egrp);

          double ell = segment_length;

          double w_avg  = (N_i[dof]/sig_t)*(1.0-exp(-sig_t*ell));
          w_avg += ((N_f[dof]-N_i[dof])/(sig_t*sig_t*ell))*
                   (1.0 - (1+sig_t*ell)*exp(-sig_t*ell));
          w_avg *= prtcl.w/ell;

          double pwl_tally_contrib = segment_length * w_avg;

          grid_tally_blocks[t].tally_local[ir]     += pwl_tally_contrib;
          grid_tally_blocks[t].tally_sqr_local[ir] += pwl_tally_contrib*
                                                      pwl_tally_contrib;
        }//for dof

        //reset for new segment
        N_i = N_f;
        prtcl.w *= exp(-sig_t*segment_length);
      }//for segment_length
    }//if make pwld
    else
      prtcl.w = w_exit;
  }//for t

}