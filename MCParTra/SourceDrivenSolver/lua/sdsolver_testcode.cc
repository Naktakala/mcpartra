#include "ChiLua/chi_lua.h"

#include "sdsolver_lua_utils.h"

namespace mcpartra
{
namespace lua_utils
{
int TestCode(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  const int num_args = lua_gettop(L);

  if (num_args != 1)
    LuaPostArgAmountError(fname, 1, num_args);

  LuaCheckNilValue(fname, L, 1);

  const int handle = lua_tointeger(L, 1);
  auto solver = mcpartra::lua_utils::GetSolverByHandle(handle, fname);

  solver->TestCode();

  return 0;
}
}//namespace lua_utils
}//namespace mcpartra