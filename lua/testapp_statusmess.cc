#include"ChiLua/chi_lua.h"
#include "chi_log.h"


int chiPrintStatus(lua_State *L)
{
  ChiLog& log = ChiLog::GetInstance();
  log.Log(LOG_0) << "Hello from here";
  return 0;
}