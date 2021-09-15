#include "sdsolver.h"

#include "chi_log.h"
#include "ChiTimer/chi_timer.h"
extern ChiTimer chi_program_timer;

//###################################################################
/**Computes PWLD transformations of the PWLD tallies.*/
void mcpartra::SourceDrivenSolver::ComputePWLDTransformations()
{
  auto& chi_log = ChiLog::GetInstance();
  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Processing PWLD Tallies.";
  if (options.make_pwld)
  {
    for (auto& cell : grid->local_cells)
    {
      auto& cell_pwl_view = pwl->GetUnitIntegrals(cell);

      MatDbl A(cell_pwl_view.GetIntV_shapeI_shapeJ());
      MatDbl Ainv = chi_math::Inverse(A);
      VecDbl b(cell_pwl_view.NumNodes(),0.0);

      auto& raw_tally = grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]];
      auto& out_tally = grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]];

      {
        if (raw_tally.empty()) continue;

        for (size_t m=0; m < num_moments; ++m)
        {
          for (size_t g=0; g < num_groups; ++g)
          {

            for (size_t i=0; i<cell_pwl_view.NumNodes(); ++i)
            {
              size_t ir = pwl->MapDOFLocal(cell, i, uk_man_pwld, m, g);
              b[i] = raw_tally.tally_global[ir]*
                     cell_pwl_view.IntV_shapeI(i);
            }//for dof

            VecDbl x = chi_math::MatMul(Ainv,b);

            for (size_t i=0; i<cell_pwl_view.NumNodes(); ++i)
            {
              size_t ir = pwl->MapDOFLocal(cell, i, uk_man_pwld, m, g);
              out_tally.tally_global[ir] = x[i];
            }//for dof

          }//for group
        }//for moment
      }//for tally

    }//for local cell lc
  }//if make_pwld

  chi_log.Log() << chi_program_timer.GetTimeString()
                << " Done processing PWLD Tallies.";
}