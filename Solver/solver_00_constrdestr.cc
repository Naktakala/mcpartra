#include "solver_montecarlon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Default constructor*/
chi_montecarlon::Solver::Solver() :
  rng0(chi_mpi.location_id)
{
  num_grps = 1;

  this->tolerance = 0.000001;

  //options
  num_particles = 1000;
  tfc_update_interval = 2000;
  mono_energy = false;
  scattering_order = 10;
  force_isotropic = false;
  group_hi_bound = -1;
  group_lo_bound = -1;
  tally_rendezvous_intvl = 100000;
  tally_multipl_factor = 1.0;
  make_pwld = false;
  uncollided_only = false;

  max_sigma = 0.0;

}

//###################################################################
/**Returns the group-mask for this solver.*/
std::pair<int,int> chi_montecarlon::Solver::GetGroupBounds()
{
  if (group_lo_bound>group_hi_bound)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::Solver: "
      << "Invalid group mask " << group_lo_bound << "->" << group_hi_bound
      << ". group_hi_bound must be >= group_lo_bound.";
    exit(EXIT_FAILURE);
  }
  int start = (group_lo_bound<0)? 0 : group_lo_bound;
  int end   = (group_hi_bound<0)? num_grps : group_lo_bound;

  return {start,end};
}