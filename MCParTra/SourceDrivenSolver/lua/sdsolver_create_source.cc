#include"ChiLua/chi_lua.h"

#include "sdsolver_lua_utils.h"

#include"SourceDrivenSolver/sdsolver.h"

#include "Sources/BoundarySource/mc_bndry_source.h"
#include "Sources/MaterialSource/mc_material_source.h"
#include "Sources/ResidualSource/mc_rmcA_source.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


namespace mcpartra
{
namespace lua_utils
{


//#############################################################################
/** Creates a simple point source at [0 0 0].
 *
\param SolverHandle int Handle to an existing montecarlo solver.
\param SourceType string Source type identifier. See SourceType below.

##_

###PropertyIndex\n
BNDRY_SRC\n
 Source on a surface boundary. Expects to be followed by the boundary number.\n\n

MATERIAL_SRC\n
 Source from material defined isotropic source. No value follows.\n\n

RESIDUAL_TYPE_A\n
 Generates source particles from a multigroup transport residual. The residual
 is computed from the flux contained in a specified field-function.
 Expects to be followed by a handle to a field-function containing flux
 moments of the approximate solution..\n\n

\return Handle int Handle to the created source.
\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonCreateSource(lua_State *L)
{
  const std::string function_name(__FUNCTION__);
  int num_args = lua_gettop(L);
  if (num_args < 2)
    LuaPostArgAmountError(__FUNCTION__, 2, num_args);

  LuaCheckNilValue(__FUNCTION__,L,1);
  LuaCheckNilValue(__FUNCTION__,L,2);

  LuaCheckNumberValue(__FUNCTION__, L, 1);
  LuaCheckStringValue(__FUNCTION__, L, 2);

  const int         solver_handle = lua_tonumber(L, 1);
  const std::string source_type   = lua_tostring(L, 2);

  auto solver = mcpartra::lua_utils::
                  GetSolverByHandle(solver_handle, function_name);

  //============================================= Boundary source
  if (source_type == "BNDRY_SRC")
  {
    if (num_args < 3)
      LuaPostArgAmountError((function_name + ": With SourceType=BOUNDARY_SOURCE"),
                            3,num_args);

    int ref_boundary = lua_tonumber(L,3);
    if (ref_boundary == 0)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid boundary number supplied in call to "
        << "chiMonteCarlonCreateSource-MC_BNDRY_SRC. Expected a positive number"
           " or a predefined identifier.";
      exit(EXIT_FAILURE);
    }


    auto new_source = new mcpartra::BoundarySource(*solver, ref_boundary);

    solver->sources.push_back(new_source);
    lua_pushinteger(L,static_cast<lua_Integer>(solver->sources.size()-1));

    chi_log.Log(LOG_0) << "MCParTra: Created boundary source.";
  }
  //============================================= Material source
  else if (source_type == "MATERIAL_SRC")
  {
    if (num_args < 2)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "MATERIAL_SRC",
                            2,num_args);

    auto new_source = new mcpartra::MaterialSource(*solver);

    solver->sources.push_back(new_source);
    lua_pushinteger(L,static_cast<lua_Integer>(solver->sources.size()-1));

    chi_log.Log(LOG_0) << "MCParTra: Created material source.";
  }
  else if (source_type == "RESIDUAL_TYPE_A")
  {
    if (num_args != 3)
      LuaPostArgAmountError("chiMonteCarlonCreateSource-"
                            "MCSrcTypes.RESIDUAL3",
                            3,num_args);

    int ff_handle = lua_tonumber(L,3);

    std::shared_ptr<chi_physics::FieldFunction> ff;
    try {
      ff = chi_physics_handler.fieldfunc_stack.at(ff_handle);
    }

    catch (std::out_of_range& o)
    {
      chi_log.Log(LOG_ALLERROR)
        << "Invalid field function handle supplied in call to "
           "chiMonteCarlonCreateSource-RESIDUAL3";
      exit(EXIT_FAILURE);
    }

    auto new_source =
      new mcpartra::ResidualSourceA(*solver, ff);

    solver->sources.push_back(new_source);
    lua_pushinteger(L,static_cast<lua_Integer>(solver->sources.size()-1));

    chi_log.Log(LOG_0) << "MCParTra: Created residual type A source.";

  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "Invalid boundary type supplied in call to chiMonteCarlonCreateSource";
    exit(EXIT_FAILURE);
  }

  return 1;
}

}//namespace lua_utils
}//namespace mcpartra