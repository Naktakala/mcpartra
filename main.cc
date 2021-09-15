#include "chi_runtime.h"

#include "ChiConsole/chi_console.h"
#include "chi_log.h"

int main(int argc, char* argv[])
{
  ChiLog&     log     = ChiLog::GetInstance();

  log.Log(LOG_0) << "MCPARTA - Monte Carlo Particle Transport";

  ChiTech::Initialize(argc,argv);

  ChiConsole& console = ChiConsole::GetInstance();

  auto L = console.consoleState;
  #include "ChiMacros/lua_register_macro.h"
  #include "SourceDrivenSolver/lua/lua_register.h"
  RegisterFunction(chiPrintStatus);

  ChiTech::RunBatch(argc,argv);

  ChiTech::Finalize();

  log.Log(LOG_0) << "MCPARTA - Execution finished";

  return 0;
}