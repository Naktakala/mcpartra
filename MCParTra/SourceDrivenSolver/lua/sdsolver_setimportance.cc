#include <ChiLua/chi_lua.h>

#include "sdsolver_lua_utils.h"

#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiPhysics/chi_physics.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"
extern ChiLog& chi_log;


namespace mcpartra
{
namespace lua_utils
{


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
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError("chiMonteCarlonSetImportances",3,num_args);

  int solver_handle = lua_tonumber(L, 1);
  int volume_hndl   = lua_tonumber(L, 2);
  double importance = lua_tonumber(L, 3);

  auto mcsolver = mcpartra::lua_utils::GetSolverByHandle(solver_handle,fname);

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

//###################################################################
/**Sets the importances of all cells from a file.

\param SolverHandle int Handle to the montecarlo solver.
\param FileName string Name of the file to read the importances from.

\author Jan
 */
int chiMonteCarlonReadImportanceMap(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError("chiMonteCarlonSetImportances",3,num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  const int solver_handle = lua_tonumber(L, 1);
  const std::string file_name = lua_tostring(L, 2);

  auto mcsolver = mcpartra::lua_utils::GetSolverByHandle(solver_handle,fname);

  mcsolver->ReadImportanceMap(file_name);

  return 0;
}

}//namespace lua_utils
}//namespace mcpartra