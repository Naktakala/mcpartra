#include <ChiLua/chi_lua.h>

#include "mcpartra_lua_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


//#############################################################################
/** Executes a MonteCarlon solver.

\param SolverHandle int Handle to the montecarlo solver.
\param PropertyIndex int Code for a specific property.



##_

###PropertyIndex\n
MC_NUM_PARTICLES\n
 Number of particles to run. Expects to be followed by an integer specifying
 the amount of particles to run. Default 1000.\n\n

MC_TFC_UPDATE_DIV\n
 Number of divisions of the number of particles to use for tally fluctuation
 chart (TFC) update. Expects to be followed by an integer specifying the number
 of bins. Default 10.\n\n

MC_MONOENERGETIC\n
 Forces the scattering out of a group to be treated like absorbtion.
 Expects to be followed by a boolean value. Default false.\n\n

MC_SCATTERING_ORDER\n
 Sets the scattering order used for building discrete scattering angles.
 Expect to be followed by an integer specifying the order. Default 10.
 Note: when this number is set greater than the scattering order available
 in the provided cross-sections then the scattering order will default to that
 available.\n\n

MC_FORCE_ISOTROPIC\n
 Flag forcing isotropic scattering. Expects to be followed by a boolean value.
 Default false.\n\n

MC_TALLY_MULTIPLICATION_FACTOR\n
 Classical global tally multiplication factor to be applied after normalization
 per source particle. Expects to be followed by a float. Default 1.0.\n\n

MC_MAKE_PWLD_SOLUTION\n
 Classical global tally multiplication factor to be applied after normalization
 per source particle. Expects to be followed by a float. Default 1.0.\n\n

MC_UNCOLLIDED_ONLY\n
 Trace source particles as uncollided particles, adjusting the weight
 instead of scattering or absorbing.\n\n

\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonSetProperty(lua_State *L)
{
  const std::string fname(__FUNCTION__);
  int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  int solver_handle = lua_tointeger(L,1);
  auto mcsolver = mcpartra::lua_utils::GetSolverByHandle(solver_handle,fname);

  //============================================= Process property index
  int property_index = lua_tonumber(L,2);
  if (property_index == mcpartra::Property::NUM_PARTICLES)
  {
    unsigned long long num_part = lua_tonumber(L,3);

    mcsolver->options.num_particles = num_part;
  }
  else if (property_index == mcpartra::Property::SCATTERING_ORDER)
  {
    int scatorder = lua_tonumber(L,3);

    mcsolver->options.scattering_order = scatorder;
  }
  else if (property_index == mcpartra::Property::MONOENERGETIC)
  {
    bool mono = lua_toboolean(L,3);

    mcsolver->options.mono_energy = mono;
  }
  else if (property_index == mcpartra::Property::FORCE_ISOTROPIC)
  {
    bool iso = lua_toboolean(L,3);

    mcsolver->options.force_isotropic = iso;
  }
  else if (property_index == mcpartra::Property::GROUP_BOUNDS)
  {
    int hi = lua_tonumber(L,3);
    int lo = lua_tonumber(L,4);

    mcsolver->options.group_hi_bound = hi;
    mcsolver->options.group_lo_bound = lo;
  }
  else if (property_index == mcpartra::Property::TALLY_MERGE_INTVL)
  {
    unsigned long long tal_merg_invtl = lua_tonumber(L,3);

    mcsolver->options.tally_rendezvous_intvl = tal_merg_invtl;
  }
  else if (property_index == mcpartra::Property::TALLY_MULTIPLICATION_FACTOR)
  {
    double tmf = lua_tonumber(L,3);

    mcsolver->options.tally_multipl_factor = tmf;
  }
  else if (property_index == mcpartra::Property::MAKE_PWLD_SOLUTION)
  {
    bool make_pwld = lua_toboolean(L,3);

    mcsolver->options.make_pwld = make_pwld;
  }
  else if (property_index == mcpartra::Property::UNCOLLIDED_ONLY)
  {
    bool unc_only = lua_toboolean(L,3);

    mcsolver->options.uncollided_only = unc_only;
  }
  else if (property_index == mcpartra::Property::NUM_UNCOLLIDED_PARTICLES)
  {
    unsigned long long num_part = lua_tonumber(L,3);

    mcsolver->options.num_uncollided_particles = num_part;
  }
  else
  {
    chi_log.Log(LOG_ALLERROR)
      << "chiMonteCarlonSetProperty:Invalid property index supplied. "
      << property_index;
    exit(EXIT_FAILURE);
  }

  return 0;
}
