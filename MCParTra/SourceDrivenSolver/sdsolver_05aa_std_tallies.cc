#include "sdsolver.h"

#include <ChiMesh/Raytrace/raytracing.h>

#include "chi_log.h"
extern ChiLog& chi_log;

typedef unsigned long long TULL;

//###################################################################
/**Makes a contribution to tallies*/
void mcpartra::SourceDrivenSolver::
  ContributeTally(const chi_mesh::Cell& cell,
                  const mcpartra::Particle &prtcl,
                  const chi_mesh::Vector3& pf)
{
  double tracklength = (pf - prtcl.pos).Norm();

  double tlw = tracklength * prtcl.w; ///< Tracklength times weight

  //============================================= Custom tallies
  for (auto& tally : custom_tallies)
  {
    if (tally.local_cell_tally_mask[cell.local_id])
    {
      for (size_t m=0; m < num_moments; ++m)
      {
        auto dof_map = uk_man_fv.MapUnknown(m,prtcl.egrp);

        double tlw_Ylm = tlw * prtcl.moment_values[m];

        tally.grid_tally.counter_local[dof_map]   += 1;
        tally.grid_tally.tally_local[dof_map]     += tlw_Ylm;
        tally.grid_tally.tally_sqr_local[dof_map] += tlw_Ylm * tlw_Ylm;
      }//for m
    }//if cell part of tally
  }//for custom tallies

  //============================================= FV Tallies
  for (unsigned int t : fv_tallies)
  {
    if (prtcl.tally_mask & (1 << t))
    {
      for (size_t m=0; m < num_moments; ++m)
      {
        int64_t dof_map = fv->MapDOFLocal(cell,/*node*/0,uk_man_fv,m,prtcl.egrp);

        double tlw_Ylm = tlw * prtcl.moment_values[m];

        grid_tally_blocks[t].counter_local[dof_map]   += 1;
        grid_tally_blocks[t].tally_local[dof_map]     += tlw_Ylm;
        grid_tally_blocks[t].tally_sqr_local[dof_map] += tlw_Ylm * tlw_Ylm;
      }//for m
    }//if tally applies
  }//for fv tallies

  //============================================= PWL Tallies
  for (unsigned int t : pwl_tallies)
  {
    if ( (prtcl.tally_mask & (1 << t)) && (options.make_pwld))
    {
      segment_lengths.clear();
      segment_lengths.push_back(tracklength);
      chi_mesh::PopulateRaySegmentLengths(*grid, cell,
                                          segment_lengths,
                                          prtcl.pos, pf,prtcl.dir);

      auto cell_pwl_view = pwl->GetCellMappingFE(cell.local_id);

      double last_segment_length = 0.0;
      for (auto segment_length : segment_lengths)
      {
        double d = last_segment_length + 0.5*segment_length;
        auto p = prtcl.pos + prtcl.dir*d;

        cell_pwl_view->ShapeValues(p, N_f);

        for (int i=0; i < cell_pwl_view->num_nodes; ++i)
        {
          for (size_t m=0; m < num_moments; ++m)
          {
            int64_t dof_map = pwl->MapDOFLocal(cell, i, uk_man_pwld,m, prtcl.egrp);

            double pwl_tlw_Ylm = segment_length * prtcl.w * N_f[i] *
                                 prtcl.moment_values[m];

            grid_tally_blocks[t].counter_local[dof_map]   += 1;
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


