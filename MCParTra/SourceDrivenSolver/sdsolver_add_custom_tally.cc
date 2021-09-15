#include "sdsolver.h"

#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/** Adds a custom volume tally to the solver.*/
size_t mcpartra::SourceDrivenSolver::
  AddCustomVolumeTally(chi_mesh::LogicalVolume &logical_volume)
{
  chi_mesh::Region*  aregion = this->regions.back();
  auto ref_grid              = aregion->GetGrid();

  if (ref_grid == nullptr)
    throw std::logic_error("Call made to " + std::string(__FUNCTION__) +
                           " without grid being assigned.");

  std::vector<bool> local_cell_tally_mask(ref_grid->local_cells.size(),false);

  for (auto& cell : ref_grid->local_cells)
    if (logical_volume.Inside(cell.centroid))
      local_cell_tally_mask[cell.local_id] = true;

  custom_tallies.emplace_back(local_cell_tally_mask);

  chi_log.Log() << "chi_montecarlon::Solver added custom tally.";

  return custom_tallies.size()-1;
}