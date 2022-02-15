#include "ChiLua/chi_lua.h"

#include "sdsolver_lua_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


namespace mcpartra
{
namespace lua_utils
{


//#############################################################################
/** Sets a property of a MonteCarlon solver.

\param SolverHandle int Handle to the montecarlo solver.
\param PropertyName string String name for a specific property.

##_

###PropertyIndex\n
"NUM_PARTICLES"\n
 Number of particles to run. Expects to be followed by an integer specifying
 the amount of particles to run. Default 1000.\n\n

"MONOENERGETIC"\n
 Forces the scattering out of a group to be treated like absorbtion.
 Expects to be followed by a boolean value. Default false.\n\n

"SCATTERING_ORDER"\n
 Sets the scattering order used for building discrete scattering angles.
 Expect to be followed by an integer specifying the order. Default 10.
 Note: when this number is set greater than the scattering order available
 in the provided cross-sections then the scattering order will default to that
 available.\n\n

"FORCE_ISOTROPIC"\n
 Flag forcing isotropic scattering. Expects to be followed by a boolean value.
 Default false.\n\n

"GROUP_BOUNDS"\n
 Bounds the simulations only between the specified groups (all inclusive).
 Expects to be followed by two integers.\n\n

"TALLY_MERGE_INTVL"\n
 Sets how many particles to run before merging tallies. Expects to be followed
 by an integer.

"TALLY_MULTIPLICATION_FACTOR"\n
 Classical global tally multiplication factor to be applied after normalization
 per source particle. Expects to be followed by a float. Default 1.0.\n\n

"MAKE_PWLD_SOLUTION"\n
 Classical global tally multiplication factor to be applied after normalization
 per source particle. Expects to be followed by a float. Default 1.0.\n\n

"UNCOLLIDED_ONLY"\n
 Trace source particles as uncollided particles, adjusting the weight
 instead of scattering or absorbing.\n\n

\ingroup LuaMonteCarlon
\author Jan*/
int chiMonteCarlonSetProperty2(lua_State *L)
{
  const std::string fname(__FUNCTION__);
  int num_args = lua_gettop(L);
  if (num_args < 3)
    LuaPostArgAmountError(__FUNCTION__, 3, num_args);

  LuaCheckNilValue(__FUNCTION__, L, 1);
  LuaCheckNilValue(__FUNCTION__, L, 2);
  LuaCheckNilValue(__FUNCTION__, L, 3);

  if (not lua_isstring(L,2))
    throw std::invalid_argument(fname + ": Argument 2 expected to be string.");

  int solver_handle = lua_tointeger(L,1);
  auto mcsolver = mcpartra::lua_utils::GetSolverByHandle(solver_handle,fname);

  //============================================= Process property index
  std::string property_index(lua_tostring(L,2));
  if (property_index == "NUM_PARTICLES")
  {
    unsigned long long num_part = lua_tonumber(L,3);

    mcsolver->options.num_particles = num_part;

    chi_log.Log() << "MCParTra: Number of particles to run set to " << num_part;
  }
  else if (property_index == "SCATTERING_ORDER")
  {
    int scatorder = lua_tonumber(L,3);

    if (scatorder > 7)
    {
      chi_log.Log() << __FUNCTION__ << ": maximum scattering order allowed "
                                       "is 7. Defaulting to 7.";
      scatorder = 7;
    }

    mcsolver->options.scattering_order = scatorder;

    chi_log.Log() << "MCParTra: Scattering order set to " << scatorder;
  }
  else if (property_index == "MONOENERGETIC")
  {
    bool mono = lua_toboolean(L,3);

    mcsolver->options.mono_energy = mono;

    chi_log.Log() << "MCParTra: Mono-energetic flag set to " << mono;
  }
  else if (property_index == "FORCE_ISOTROPIC")
  {
    bool iso = lua_toboolean(L,3);

    mcsolver->options.force_isotropic = iso;

    chi_log.Log() << "MCParTra: Force-isotropic flag set to " << iso;
  }
  else if (property_index == "GROUP_BOUNDS")
  {
    int hi = lua_tonumber(L,3);
    int lo = lua_tonumber(L,4);

    mcsolver->options.group_hi_bound = hi;
    mcsolver->options.group_lo_bound = lo;

    chi_log.Log() << "MCParTra: Group-bounds set from " << lo << " to " << hi;
  }
  else if (property_index == "TALLY_MERGE_INTVL")
  {
    unsigned long long tal_merg_invtl = lua_tonumber(L,3);

    mcsolver->options.tally_rendezvous_intvl = tal_merg_invtl;

    chi_log.Log() << "MCParTra: Tally merge interval set to " << tal_merg_invtl;
  }
  else if (property_index == "TALLY_MULTIPLICATION_FACTOR")
  {
    double tmf = lua_tonumber(L,3);

    mcsolver->options.tally_multipl_factor = tmf;

    chi_log.Log() << "MCParTra: Tally multiplication factor set to " << tmf;
  }
  else if (property_index == "MAKE_PWLD_SOLUTION")
  {
    bool make_pwld = lua_toboolean(L,3);

    mcsolver->options.make_pwld = make_pwld;

    chi_log.Log() << "MCParTra: Flag to make PWLD-solution set to " << make_pwld;
  }
  else if (property_index == "UNCOLLIDED_ONLY")
  {
    bool unc_only = lua_toboolean(L,3);

    mcsolver->options.uncollided_only = unc_only;

    chi_log.Log() << "MCParTra: Flag to run uncollided only, set to " << unc_only;
  }
  else if (property_index == "RUN_TAPE_BASE_NAME")
  {
    if (not lua_isstring(L,3))
      throw std::invalid_argument(fname + ": Argument 3 expected to be a string.");

    std::string base_name(lua_tostring(L,3));

    mcsolver->options.run_tape_base_name = base_name;
    mcsolver->options.write_run_tape = true;

    chi_log.Log() << "MCParTra: Run-tape base name set to " << base_name
                  << ", and flag to write run tape, set to true.";
  }
  else if (property_index == "IMPORTANCE_DURING_RAYTRACING")
  {
    LuaCheckBoolValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    mcsolver->options.importances_during_raytracing = flag;
    if (flag) mcsolver->options.apply_source_importance_sampling = false;

    chi_log.Log() << "MCParTra: Flag to use importances during raytracing set "
                  << "to " << flag << ".";
  }
  else if (property_index == "APPLY_SOURCE_IMPORTANCE_SAMPLING")
  {
    LuaCheckBoolValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    mcsolver->options.apply_source_importance_sampling = flag;
    if (flag) mcsolver->options.importances_during_raytracing = false;

    chi_log.Log() << "MCParTra: Flag to apply importance biased source "
                  << "sampling set to " << flag << ".";
  }
  else if (property_index == "APPLY_SOURCE_ANGULAR_BIASING")
  {
    LuaCheckBoolValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    mcsolver->options.apply_source_angular_biasing = flag;
    if (flag) mcsolver->options.importances_during_raytracing = false;

    chi_log.Log() << "MCParTra: Flag to apply angular importance biased source "
                  << "sampling set to " << flag << ".";
  }
  else if (property_index == "PRINT_TFCS")
  {
    LuaCheckBoolValue(fname, L, 3);

    const bool flag = lua_toboolean(L, 3);

    mcsolver->options.print_TFC = flag;

    chi_log.Log() << "MCParTra: Flag to print tally fluctuation charts set"
                     " to " << flag << ".";
  }
  else if (property_index == "RESIDUAL_SRC_FF_OPTION")
  {
    LuaCheckNumberValue(fname, L, 3);

    unsigned int option = lua_tointeger(L, 3);

    using RFOPT = SourceDrivenSolver::ResidSrcFFOption;

    auto option_set = RFOPT::DISCONTINUOUS_Q1;

    if      (option == 0) option_set = RFOPT::DISCONTINUOUS_Q1;
    else if (option == 1) option_set = RFOPT::DISCONTINUOUS_Q0;
    else if (option == 2) option_set = RFOPT::CONTINUOUS_Q1;
    else
      throw std::invalid_argument(fname + ": Invalid value for option "
                                          "RESIDUAL_SRC_FF_OPTION.");

    mcsolver->options.resid_src_ff_option = option_set;

    chi_log.Log() << "MCParTra: Residual approximate solution treatment set to "
                  << option << ".";
  }
  else if (property_index == "RESIDUAL_SRC_NY")
  {
    LuaCheckNumberValue(fname, L, 3);
    const int N_y = lua_tonumber(L, 3);

    mcsolver->options.resid_src_integration_N_y = N_y;

    chi_log.Log() << "MCParTra: Residual-source number of integration points "
                     "set to " << N_y << ".";
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

}//namespace lua_utils
}//namespace mcpartra
