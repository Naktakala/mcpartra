#include"ChiLua/chi_lua.h"

#include "SourceDrivenSolver/sdsolver.h"
#include "ChiMesh/MeshHandler/chi_meshhandler.h"
#include "ChiMesh/LogicalVolume/chi_mesh_logicalvolume.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


namespace mcpartra
{
namespace lua_utils
{


//#############################################################################
/** Designates a custom volumetric tally bounded by a logical volume.
 *
\param solver_handle int Handle to the solver.
\param lv_handle     int Handle to the logical volume.

\return Handle to the newly created custom tally.*/
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

  auto mc_solver = dynamic_cast<mcpartra::SourceDrivenSolver*>(solver);
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

//#############################################################################
/** Retrieves the value of a custom tally.
 *
\param solver_handle int Handle to the solver.
\param tally_handle  int Handle to the custom tally on the solver.
\param moment        int The specific moment required.
\param group         int The specific group required.

\return Returns the tally value and its 1-sigma stddev.*/
int chiMonteCarlonGetCustomVolumeTallyValue(lua_State *L)
{
  //============================================= Check arguments
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 4)
    LuaPostArgAmountError(fname,4,num_args);

  LuaCheckNilValue(fname,L,1);
  LuaCheckNilValue(fname,L,2);
  LuaCheckNilValue(fname,L,3);
  LuaCheckNilValue(fname,L,4);

  //============================================= Parse arguments
  //----------------------------------- Solver
  const int solver_index = lua_tonumber(L,1);

  chi_physics::Solver* solver;
  try{
    solver = chi_physics_handler.solver_stack.at(solver_index);
  }
  catch(const std::out_of_range& o){
    throw std::invalid_argument(fname + ": Invalid solver handle in call to.");
  }

  const auto mc_solver = dynamic_cast<mcpartra::SourceDrivenSolver*>(solver);
  if (not mc_solver)
    throw std::invalid_argument(fname + ": Solver pointed to does not have"
                                " type: mcpartra::SourceDrivenSolver.");

  const auto uk_man_fv = mc_solver->GetUnknownManagerFV();

  //----------------------------------- Logical Surface
  const int tally_handle = lua_tointeger(L, 2);
  const int moment       = lua_tointeger(L, 3);
  const int group        = lua_tointeger(L, 4);

  const auto& tally = mc_solver->GetCustomVolumeTally(tally_handle);

  if (tally.tally_fluctation_chart.empty())
    throw std::logic_error(fname + ": Custom tally has no populated values. "
                                   "This might be because the simulation has "
                                   "not been executed.");

  const auto tfc_entry = tally.tally_fluctation_chart.back();

  const unsigned int dof_map = uk_man_fv.MapUnknown(moment, group);

  double tally_value = 0.0;
  double tally_sigma = 1.0;
  try {
    tally_value = tfc_entry.average.at(dof_map);
    tally_sigma = tfc_entry.sigma.at(dof_map);
  }
  catch (const std::out_of_range& oor)
  {
    throw std::logic_error(fname + ": Custom tally has been accessed at a "
                                   "DOF-index outside that which is available. "
                                   "Check the number of available groups and"
                                   " moments.");
  }

  lua_pushnumber(L, tally_value);
  lua_pushnumber(L, tally_sigma);
  return 2;
}

}//namespace lua_utils
}//namespace mcpartra