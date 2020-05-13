#include "solver_montecarlon.h"

#include <chi_mpi.h>

extern ChiMPI& chi_mpi;

//###################################################################
/**Default constructor*/
chi_montecarlon::Solver::Solver() :
  rng0(chi_mpi.location_id,0)
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

  max_relative_error = 0.0;

}