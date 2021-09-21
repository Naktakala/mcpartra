#include "sdsolver_lua_utils.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

mcpartra::SourceDrivenSolver* mcpartra::lua_utils::
  GetSolverByHandle(int handle, const std::string& function_name)
{
  mcpartra::SourceDrivenSolver* mcpartra_solver;
  try{

    mcpartra_solver = dynamic_cast<mcpartra::SourceDrivenSolver*>(
      chi_physics_handler.solver_stack.at(handle));

    if (not mcpartra_solver)
      throw std::logic_error(function_name + ": Invalid solver at given handle (" +
                             std::to_string(handle) + "). "
                             "The solver is not of type mcpartra::Solver.");
  }//try
  catch(const std::out_of_range& o) {
    throw std::logic_error(function_name + ": Invalid solver-handle (" +
                           std::to_string(handle) + ").");
  }

  return mcpartra_solver;
}

void mcpartra::lua_utils::RegisterLuaEntities(lua_State *L)
{
  lua_register(L, "chiMonteCarlonCreateSolver",
               mcpartra::lua_utils::chiMonteCarlonCreateSolver);
  lua_register(L, "chiMonteCarlonCreateSource",
               mcpartra::lua_utils::chiMonteCarlonCreateSource);

  lua_register(L, "chiMonteCarlonReadRuntape",
               mcpartra::lua_utils::chiMonteCarlonReadRuntape);
  lua_register(L, "chiMonteCarlonWriteLBSFluxMoments",
               mcpartra::lua_utils::chiMonteCarlonWriteLBSFluxMoments);
  lua_register(L, "chiMonteCarlonSetImportances",
               mcpartra::lua_utils::chiMonteCarlonSetImportances);
  lua_register(L, "chiMonteCarlonSetProperty2",
               mcpartra::lua_utils::chiMonteCarlonSetProperty2);
  lua_register(L, "chiMonteCarlonAddCustomVolumeTally",
               mcpartra::lua_utils::chiMonteCarlonAddCustomVolumeTally);
}