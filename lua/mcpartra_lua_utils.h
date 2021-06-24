#ifndef MCPARTRA_LUA_UTILS_H
#define MCPARTRA_LUA_UTILS_H

#include "../Solver/solver_montecarlon.h"

namespace mcpartra
{

namespace lua_utils
{
  mcpartra::Solver*
    GetSolverByHandle(int handle,const std::string& function_name);

}//namespace lua_utils

}//namespace mcpartra

#endif //MCPARTRA_LUA_UTILS_H