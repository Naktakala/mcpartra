#include "ChiLua/chi_lua.h"

#include "mcpartra_lua_utils.h"

#include "chi_log.h"

//###################################################################
/**Reads a runtape file given by the file_name.
 *
\param handle int Handle to the solver.
\param file_name string Complete path to the runtape file.

*/
int chiMonteCarlonReadRuntape(lua_State* L)
{
  const std::string fname = __FUNCTION__;
  int num_args = lua_gettop(L);
  if (num_args != 2)
    LuaPostArgAmountError(fname.c_str(), 2, num_args);

  LuaCheckNilValue(fname.c_str(), L, 1);
  LuaCheckNilValue(fname.c_str(), L, 2);

  if (not lua_isstring(L, 2))
    throw std::invalid_argument(fname + ": Argument 2 expected to be string.");

  int solver_handle = lua_tonumber(L, 1);
  auto mcsolver = mcpartra::lua_utils::GetSolverByHandle(solver_handle, fname);

  std::string file_name = lua_tostring(L, 2);

  mcsolver->ReadRunTape(file_name);

  return 0;
}