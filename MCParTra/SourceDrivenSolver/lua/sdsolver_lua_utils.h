#ifndef MCPARTRA_LUA_UTILS_H
#define MCPARTRA_LUA_UTILS_H

#include "SourceDrivenSolver/sdsolver.h"

namespace mcpartra
{

namespace lua_utils
{
  mcpartra::SourceDrivenSolver*
    GetSolverByHandle(int handle,const std::string& function_name);

  int chiMonteCarlonCreateSolver(lua_State* L);
  int chiMonteCarlonCreateSource(lua_State* L);

  int chiMonteCarlonReadRuntape(lua_State* L);
  int chiMonteCarlonWriteLBSFluxMoments(lua_State* L);
  int chiMonteCarlonSetImportances(lua_State* L);
  int chiMonteCarlonReadImportanceMap(lua_State* L);
  int chiMonteCarlonExportImportanceMap(lua_State* L);
  int chiMonteCarlonSetProperty2(lua_State* L);
  int chiMonteCarlonAddCustomVolumeTally(lua_State* L);

  void RegisterLuaEntities(lua_State* L);
}//namespace lua_utils

}//namespace mcpartra

#endif //MCPARTRA_LUA_UTILS_H