#ifndef MCPARTRA_LUA_UTILS_H
#define MCPARTRA_LUA_UTILS_H

#include "SourceDrivenSolver/sdsolver.h"

namespace mcpartra
{

namespace lua_utils
{
  mcpartra::SourceDrivenSolver*
    GetSolverByHandle(int handle,const std::string& function_name);

}//namespace lua_utils

}//namespace mcpartra

#endif //MCPARTRA_LUA_UTILS_H