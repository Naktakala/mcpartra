#include <ChiLua/chi_lua.h>

#include"SourceDrivenSolver/sdsolver.h"

#include <ChiPhysics/chi_physics.h>

extern ChiPhysics&  chi_physics_handler;

#include "chi_log.h"

extern ChiLog& chi_log;

//#############################################################################
/** Executes a MonteCarlon solver.

\param SolverHandle int Handle to the montecarlo solver

\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonExecute(lua_State *L)
{
  int num_args = lua_gettop(L);
  if (num_args != 1)
    LuaPostArgAmountError("chiMonteCarlonExecute",1,num_args);

  chi_physics::Solver* solver = nullptr;
  try{
    solver = chi_physics_handler.solver_stack.at(lua_tonumber(L,1));
  }
  catch (const std::out_of_range& o)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonExecute: Invalid solver handle. "
      << lua_tonumber(L,1);
    exit(EXIT_FAILURE);
  }

  if (typeid(*solver) == typeid(mcpartra::SourceDrivenSolver))
  {
    auto mcsolver = (mcpartra::SourceDrivenSolver*)solver;
    mcsolver->Execute();
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonExecute: Solver pointed to by solver handle is "
      << " not a MonteCarlo solver.";
    exit(EXIT_FAILURE);
  }
  return 0;
}