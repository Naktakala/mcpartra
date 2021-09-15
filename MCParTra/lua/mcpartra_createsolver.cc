#include"ChiLua/chi_lua.h"

#include"SourceDrivenSolver/sdsolver.h"

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

#include"ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


//#############################################################################
/** Creates a MonteCarlon solver.

\return Handle int Handle to the created solver.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSolver(lua_State *L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);

  if (num_args >= 1) LuaCheckNilValue(fname, L, 1);
  if (num_args == 2) LuaCheckNilValue(fname, L, 2);

  int                seed = chi_mpi.location_id;
  std::string solver_name = "MCParTra";

  if (num_args >= 1)
  {
    LuaCheckNumberValue(fname, L, 1);
    seed = lua_tointeger(L, 1);
  }

  if (num_args == 2)
  {
    LuaCheckStringValue(fname, L, 2);
    solver_name = lua_tostring(L,2);
  }

  auto new_solver = new mcpartra::SourceDrivenSolver(seed, solver_name);

  chi_physics_handler.solver_stack.push_back(new_solver);

  auto stack_size = chi_physics_handler.solver_stack.size();
  lua_pushinteger(L,static_cast<lua_Integer>(stack_size)-1);

  return 1;
}
