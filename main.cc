#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

#include "MCParTra/SourceDrivenSolver/lua/sdsolver_lua_utils.h"

int main(int argc, char* argv[])
{
  ChiLog&     log     = ChiLog::GetInstance();

  log.Log(LOG_0) << "MCPARTA - Monte Carlo Particle Transport";

  ChiTech::Initialize(argc,argv);

  auto& lua_console = ChiConsole::GetInstance();

  mcpartra::lua_utils::RegisterLuaEntities(lua_console.consoleState);

  ChiTech::RunBatch(argc,argv);

  ChiTech::Finalize();

  log.Log(LOG_0) << "MCPARTA - Execution finished";

  return 0;
}