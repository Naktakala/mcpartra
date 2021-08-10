#include "solver_montecarlon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Default constructor*/
mcpartra::Solver::Solver(const std::string& text_name) :
  chi_physics::Solver(text_name),
  rng0(chi_mpi.location_id)
{
  chi_log.Log() << "MCParTra: Solver created with seed " << chi_mpi.location_id;
}

//###################################################################
/**Constructor with seed.*/
mcpartra::Solver::Solver(int seed, const std::string& text_name) :
  chi_physics::Solver(text_name),
  rng0(seed)
{
  chi_log.Log() << "MCParTra: Solver created with seed " << seed;
}