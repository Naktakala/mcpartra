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
  uint64_t cell_local_ind = cell.local_id;

  double tracklength = (pf - prtcl.pos).Norm();

  double tlw = tracklength * prtcl.w; ///< Tracklength times weight

  if (std::isnan(tracklength))
  {
    chi_log.Log(LOG_ALLERROR)
      << "Tracklength corruption."
      << " pos  " << prtcl.pos.PrintS()
      << " posf " << pf.PrintS();
    exit(EXIT_FAILURE);
  }

  //============================================= Custom tallies

  for (auto& tally : custom_tallies)
  {
    if (tally.local_cell_tally_mask[cell.local_id])
    {
      for (int m=0; m<num_moms; ++m)
      {
        auto dof_map = uk_man_fv.MapUnknown(m,prtcl.egrp);

        double tlw_Ylm = tlw * prtcl.moment_values[m];

        tally.grid_tally.tally_local[dof_map]     += tlw_Ylm;
        tally.grid_tally.tally_sqr_local[dof_map] += tlw_Ylm * tlw_Ylm;
      }//for m
    }//if cell part of tally
  }//for custom tallies

  //============================================= FV Tallies
  for (int t : fv_tallies)
  {
    if (prtcl.tally_mask & (1 << t))
    {
      for (int m=0; m<num_moms; ++m)
      {
        int64_t dof_map = fv->MapDOFLocal(cell,/*node*/0,uk_man_fv,m,prtcl.egrp);

        double tlw_Ylm = tlw * prtcl.moment_values[m];

        grid_tally_blocks[t].tally_local[dof_map]     += tlw_Ylm;
        grid_tally_blocks[t].tally_sqr_local[dof_map] += tlw_Ylm * tlw_Ylm;
      }//for m
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

        for (int i=0; i < cell_pwl_view->num_nodes; ++i)
        {
          for (int m=0; m<num_moms; ++m)
          {
            int64_t dof_map = pwl->MapDOFLocal(cell, i, uk_man_pwld,/*m*/0, prtcl.egrp);
            double pwl_tlw_Ylm = segment_length * prtcl.w * N_f[i];

            grid_tally_blocks[t].tally_local[dof_map]     += pwl_tlw_Ylm;
            grid_tally_blocks[t].tally_sqr_local[dof_map] += pwl_tlw_Ylm *
                                                             pwl_tlw_Ylm;
          }//for m
        }//for node

        last_segment_length += segment_length;
      }//for segment_length
    }//if tally applies
  }//for t

}



