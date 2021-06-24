#include <ChiLua/chi_lua.h>

#include"../Solver/solver_montecarlon.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

extern ChiPhysics&  chi_physics_handler;

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**Sets the importances of certain cells using a logical volume.

\param SolverHandle int Handle to the montecarlo solver.
\param LogicalVolumeHandle int Handle to the logical volume to use for this
                               operation.
\param Importance float Value of the importance to set to all cells
                  within the volume.


\author Jan
 */
int chiMonteCarlonSetImportances(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError("chiMonteCarlonSetImportances",3,num_args);

  chi_physics::Solver* solver = nullptr;
  try{
    solver = chi_physics_handler.solver_stack.at(lua_tonumber(L,1));
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetImportances: Invalid solver handle. "
      << lua_tonumber(L,1);
    exit(EXIT_FAILURE);
  }

  chi_montecarlon::Solver* mcsolver;
  if (typeid(*solver) == typeid(chi_montecarlon::Solver))
    mcsolver = (chi_montecarlon::Solver*)solver;
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetProperty: Solver pointed to by solver handle is "
      << " not a MonteCarlo solver.";
    exit(EXIT_FAILURE);
  }

  int volume_hndl = lua_tonumber(L,2);
  double importance = lua_tonumber(L,3);

  //============================================= Get current mesh handler
  chi_mesh::MeshHandler* cur_hndlr = chi_mesh::GetCurrentHandler();

  if (volume_hndl >= cur_hndlr->logicvolume_stack.size())
  {
    chi_log.Log(LOG_ALLERROR) << "Invalid logical volume specified in "
                                 "chiVolumeMesherSetProperty("
                                 "MATID_FROMLOGICAL...";
    exit(EXIT_FAILURE);
  }

  chi_mesh::LogicalVolume* volume_ptr =
    cur_hndlr->logicvolume_stack[volume_hndl];

  //============================================= Get the grid
  auto grid = cur_hndlr->GetGrid();

  if (mcsolver->local_cell_importance_setting.empty())
  {
    mcsolver->local_cell_importance_setting.clear();
    mcsolver->local_cell_importance_setting.resize(grid->local_cells.size(),1.0);
  }

  for (auto& cell : grid->local_cells)
    if (volume_ptr->Inside(cell.centroid))
      mcsolver->local_cell_importance_setting[cell.local_id] = importance;

  return 0;
}