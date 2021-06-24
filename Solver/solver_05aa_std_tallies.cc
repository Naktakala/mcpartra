#include "solver_montecarlon.h"

#include <ChiMesh/Raytrace/raytracing.h>

#include <chi_log.h>
extern ChiLog& chi_log;

typedef unsigned long long TULL;

//###################################################################
/**Makes a contribution to tallies*/
void mcpartra::Solver::
  ContributeTally(mcpartra::Particle &prtcl, const chi_mesh::Vector3& pf)
{
  auto& cell = grid->local_cells[prtcl.cur_cell_local_id];
  int cell_local_ind = cell.local_id;

  int ir_cell = fv->MapDOFLocal(cell, 0, uk_man_fv,/*m*/0, prtcl.egrp);

  double tracklength = (pf - prtcl.pos).Norm();

  double tally_contrib = tracklength*prtcl.w;

  if (std::isnan(tracklength))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Tracklength corruption."
      << " pos  " << prtcl.pos.PrintS()
      << " posf " << pf.PrintS();
    exit(EXIT_FAILURE);
  }
  //============================================= Custom tallies
  auto ir_ctally = uk_man_fv.MapUnknown(/*m*/0,/*g*/0);
  for (auto& tally : custom_tallies)
  {
    if (tally.local_cell_tally_mask[cell.local_id])
    {
      tally.grid_tally.tally_local[ir_ctally]     += tally_contrib;
      tally.grid_tally.tally_sqr_local[ir_ctally] += tally_contrib*
                                                     tally_contrib;
    }
  }

  //============================================= FV Tallies
  for (int t : fv_tallies)
  {
    if (prtcl.tally_mask & (1 << t))
    {
      grid_tally_blocks[t].tally_local[ir_cell]     += tally_contrib;
      grid_tally_blocks[t].tally_sqr_local[ir_cell] += tally_contrib*
                                                       tally_contrib;
    }//if tally applies
  }//for fv tallies

  //============================================= PWL Tallies
  for (int t : pwl_tallies)
  {
    if ( (prtcl.tally_mask & (1 << t)) && (options.make_pwld))
    {
      segment_lengths.clear();
      segment_lengths.push_back(tracklength);
      chi_mesh::PopulateRaySegmentLengths(*grid, cell,
                                          segment_lengths,
                                          prtcl.pos, pf,prtcl.dir);

      auto cell_pwl_view = pwl->GetCellMappingFE(cell_local_ind);

      double last_segment_length = 0.0;
      for (auto segment_length : segment_lengths)
      {
        double d = last_segment_length + 0.5*segment_length;
        auto p = prtcl.pos + prtcl.dir*d;

        cell_pwl_view->ShapeValues(p, N_f);

        for (int dof=0; dof<cell_pwl_view->num_nodes; dof++)
        {
          int ir = pwl->MapDOFLocal(cell, dof, uk_man_pwld,/*m*/0, prtcl.egrp);
          double pwl_tally_contrib = segment_length * prtcl.w * N_f[dof];

          grid_tally_blocks[t].tally_local[ir]     += pwl_tally_contrib;
          grid_tally_blocks[t].tally_sqr_local[ir] += pwl_tally_contrib*
                                                      pwl_tally_contrib;
        }//for dof

        last_segment_length += segment_length;
      }//for segment_length
    }//if tally applies
  }//for t

}



