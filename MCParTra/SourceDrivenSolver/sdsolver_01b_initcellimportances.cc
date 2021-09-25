#include "sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Set cell importances.*/
void mcpartra::SourceDrivenSolver::InitCellImportances()
{
  chi_log.Log() << "MCParTra: Initializing cell importances";

  size_t num_local_cells = grid->local_cells.size();
  local_cell_importance.assign(num_local_cells,1.0);

  if (not local_cell_importance_setting.empty())
    local_cell_importance = local_cell_importance_setting;
}