#include "mcpartra_lua_utils.h"

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


//  chi_physics::Solver* generic_solver;
//  try {
//    generic_solver = chi_physics_handler.solver_stack.at(handle);
//  }
//  catch(const std::out_of_range& o) {
//    throw std::logic_error(function_name + ": Invalid solver-handle (" +
//                           std::to_string(handle) + ").");
//  }
//
//  auto mcpartra_solver = dynamic_cast<mcpartra::Solver*>(generic_solver);
//
//  if (not mcpartra_solver)
//    throw std::logic_error(function_name + ": Invalid solver at given handle (" +
//                                           std::to_string(handle) + "). "
//                                           "The solver is not of type "
//                                           "mcpartra::Solver. Type="+
//                             typeid(*generic_solver).name());

  return mcpartra_solver;
}