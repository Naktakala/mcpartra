#include "ChiLua/chi_lua.h"

#include "sdsolver_lua_utils.h"

namespace mcpartra
{
namespace lua_utils
{


//###################################################################
/**Writes LBS formatted flux moments from the pwl tally.*/
int chiMonteCarlonWriteLBSFluxMoments(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname, 2, num_args);

  LuaCheckNilValue(fname, L, 1);
  LuaCheckNilValue(fname, L, 2);

  LuaCheckIntegerValue(fname, L, 1);
  LuaCheckStringValue(fname, L, 2);

  const int solver_handle = lua_tointeger(L, 1);
  const std::string file_name = lua_tostring(L, 2);

  auto solver = mcpartra::lua_utils::GetSolverByHandle(solver_handle,fname);

  solver->WriteLBSFluxMoments(file_name);

  return 0;
}

}//namespace lua_utils
}//namespace mcpartra