#include "sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Default constructor*/
mcpartra::SourceDrivenSolver::SourceDrivenSolver(const std::string& text_name) :
  chi_physics::Solver(text_name,{{"num_particles",          int64_t(1000)},
                                 {"tally_rendezvous_intvl", int64_t(100000)},
                                 {"mono_energy",            false},
                                 {"scattering_order",       int64_t(1)},
                                 {"force_isotropic",        false},
                                 {"group_hi_bound",         int64_t(-1)},
                                 {"group_lo_bound",         int64_t(-1)},
                                 {"tally_multipl_factor",   1.0},
                                 {"make_pwld",              false},
                                 {"uncollided_only",        false},
                                 {"write_run_tape",         false},
                                 {"run_tape_base_name",     ""}}),
  rng0(chi_mpi.location_id)
{
  chi_log.Log() << "MCParTra: Solver created with seed " << chi_mpi.location_id;
}

//###################################################################
/**Constructor with seed.*/
mcpartra::SourceDrivenSolver::SourceDrivenSolver(int seed, const std::string& text_name) :
chi_physics::Solver(text_name,{{"num_particles",          int64_t(1000)},
                               {"tally_rendezvous_intvl", int64_t(100000)},
                               {"mono_energy",            false},
                               {"scattering_order",       int64_t(1)},
                               {"force_isotropic",        false},
                               {"group_hi_bound",         int64_t(-1)},
                               {"group_lo_bound",         int64_t(-1)},
                               {"tally_multipl_factor",   1.0},
                               {"make_pwld",              false},
                               {"uncollided_only",        false},
                               {"write_run_tape",         false},
                               {"run_tape_base_name",     ""}}),
  rng0(seed)
{
  chi_log.Log() << "MCParTra: Solver created with seed " << seed;
}