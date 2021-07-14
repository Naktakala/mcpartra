#include"ChiLua/chi_lua.h"

#include"../Solver/solver_montecarlon.h"

#include"ChiPhysics/chi_physics.h"

extern ChiPhysics&  chi_physics_handler;


//#############################################################################
/** Creates a MonteCarlon solver.

\return Handle int Handle to the created solver.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSolver(lua_State *L)
{
  int num_args = lua_gettop(L);

  if (num_args == 0)
  {
    auto new_solver = new mcpartra::Solver;

    chi_physics_handler.solver_stack.push_back(new_solver);

    auto stack_size = chi_physics_handler.solver_stack.size();
    lua_pushinteger(L,static_cast<lua_Integer>(stack_size)-1);

    return 1;
  }
  else
  {
    int seed = lua_tointeger(L,1);
    auto new_solver = new mcpartra::Solver(seed);

    chi_physics_handler.solver_stack.push_back(new_solver);

    auto stack_size = chi_physics_handler.solver_stack.size();
    lua_pushinteger(L,static_cast<lua_Integer>(stack_size)-1);

    return 1;
  }
}
