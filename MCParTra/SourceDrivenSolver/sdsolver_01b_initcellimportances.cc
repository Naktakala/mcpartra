#include "sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Set cell importances.*/
void mcpartra::SourceDrivenSolver::InitCellImportances()
{
  chi_log.Log() << "MCParTra: Initializing cell importances";

  size_t num_local_cells = grid->local_cells.size();
  if (local_cell_importance_setting.empty())
  {
    local_cell_importance.clear();
    local_cell_importance.resize(num_local_cells,1.0);
  }
  else
    local_cell_importance = local_cell_importance_setting;
}