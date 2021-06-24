#include"ChiLua/chi_lua.h"

#include "../Solver/solver_montecarlon.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//#############################################################################
/** Creates a simple point source at [0 0 0].*/
int chiMonteCarlonAddCustomVolumeTally(lua_State *L)
{
  //============================================= Check arguments
  int num_args = lua_gettop(L);

  if (num_args != 2)
    LuaPostArgAmountError(__FUNCTION__,2,num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  //============================================= Parse arguments
  //----------------------------------- Solver
  int solver_index = lua_tonumber(L,1);

  chi_physics::Solver* solver;
  try{
    solver = chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(const std::out_of_range& o){
    throw std::invalid_argument("Invalid solver handle in call to " +
                                std::string(__FUNCTION__) + ".");
  }

  auto mc_solver = dynamic_cast<chi_montecarlon::Solver*>(solver);
  if (not mc_solver)
    throw std::invalid_argument("Solver pointed to in call to " +
                                std::string(__FUNCTION__) + " does not have"
                                " type: chi_montecarlon::Solver.");

  //----------------------------------- Logical Surface
  int logical_volume_handle = lua_tonumber(L, 2);

  auto handler = chi_mesh::GetCurrentHandler();

  chi_mesh::LogicalVolume* logical_vol;
  try{
    logical_vol = handler->logicvolume_stack.at(logical_volume_handle);
  }
  catch(const std::out_of_range& o){
    throw std::invalid_argument("Invalid logical volume handle in call to " +
                                std::string(__FUNCTION__) + ".");
  }

  auto new_handle = mc_solver->AddCustomVolumeTally(*logical_vol);

  lua_pushnumber(L,new_handle);
  return 1;
}