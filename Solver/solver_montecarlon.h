#ifndef _solver_montecarlon
#define _solver_montecarlon

#include"../chi_montecarlon.h"
#include "../chi_montecarlon_particle.h"
#include"ChiPhysics/SolverBase/chi_solver.h"
#include "../RandomNumberGenerator/montecarlon_rng.h"
#include "../Source/mc_base_source.h"
#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <FiniteVolume/fv.h>
#include <PiecewiseLinear/pwl.h>
#include <ChiMath/chi_math.h>

#include <ChiMesh/Raytrace/raytracing.h>

#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>

namespace chi_montecarlon
{
  enum Property{
    NUM_PARTICLES               = 1,
    TFC_UPDATE_INTVL            = 2,
    MONOENERGETIC               = 3,
    SCATTERING_ORDER            = 4,
    FORCE_ISOTROPIC             = 5,
    GROUP_BOUNDS                = 6,
    TALLY_MERGE_INTVL           = 7,
    TALLY_MULTIPLICATION_FACTOR = 8,
    MAKE_PWLD_SOLUTION          = 9,
    UNCOLLIDED_ONLY             = 10,
    NUM_UNCOLLIDED_PARTICLES    = 11
  };
}

//######################################################### Class def
/**Monte Carlo neutron particle solver.*/
class chi_montecarlon::Solver : public chi_physics::Solver
{
private:
  chi_mesh::MeshContinuum*              grid;
public:
  SpatialDiscretization_FV*             fv_discretization;
  SpatialDiscretization_PWL*            pwl_discretization;
private:
  std::vector<int>                      matid_xs_map;
  std::vector<int>                      matid_q_map;

  std::vector<unsigned long long>       batch_sizes;
  std::vector<unsigned long long>       batch_sizes_per_loc;

  std::vector<unsigned long long>       uncollided_batch_sizes;
  std::vector<unsigned long long>       uncollided_batch_sizes_per_loc;
public:
  int                                   num_grps;
private:
  //FV tallies
  std::vector<double>                   phi_tally_contrib;
  std::vector<double>                   phi_tally;
  std::vector<double>                   phi_tally_sqr;

public:
  std::vector<double>                   phi_global_initial_value;
private:
  std::vector<double>                   phi_global;
  std::vector<double>                   phi_global_tally_sqr;

  std::vector<double>                   phi_local_relsigma;

  std::vector<double>                   phi_uncollided_rmc;

  //PWL tallies
  int                                   num_moms;
  std::vector<double>                   phi_pwl_tally_contrib;
  std::vector<double>                   phi_pwl_tally;
  std::vector<double>                   phi_pwl_tally_sqr;

  std::vector<double>                   phi_pwl_global;
  std::vector<double>                   phi_pwl_global_tally_sqr;

  std::vector<double>                   phi_pwl_local_relsigma;
public:
  std::vector<double>                   phi_pwl_uncollided_rmc;
private:
  std::vector<int>                      local_cell_pwl_dof_array_address;

  std::map<int,int>                     cell_neighbor_nonlocal_local_id;

  //======================== Variance reduction quantities
public:
  std::vector<double>                   local_cell_importance_setting;
private:
  std::vector<double>                   local_cell_importance;

  //======================== RMC quantities
public:
  std::vector<double>                   cell_residual_cdf;
  double                                domain_volume = 0.0;
private:
  std::vector<double> segment_lengths;
  std::vector<double> N_f,N_i;
  std::vector<chi_mesh::Vector3> Grad;
  //runtime quantities
  size_t                                current_batch;
  unsigned long long                    nps;
  unsigned long long                    nps_global;
  double                                max_relative_error;
  double                                max_relative_error2;
  double                                max_relative_error3;
  chi_math::CDFSampler*                 surface_source_sampler;

  MPI_Datatype                          mpi_prtcl_data_type;
  std::vector<Particle>                 outbound_particle_bank;
  std::vector<Particle>                 inbound_particle_bank;
  std::vector<Particle>                 particle_source_bank;
  unsigned int                          total_outbound_bank_size=0;
  unsigned int                          total_inbound_bank_size=0;

public:
  RandomNumberGenerator                 rng0;
  std::vector<chi_montecarlon::Source*> sources;

  //options
  double                                tolerance = 0.000001;
  unsigned long long                    num_particles = 1000;
  int                                   tfc_update_interval = 2000;
  bool                                  mono_energy = false;
  int                                   scattering_order = 10;
  bool                                  force_isotropic = false;
  int                                   group_hi_bound = -1;
  int                                   group_lo_bound = -1;
  unsigned long long                    tally_rendezvous_intvl = 100000;
  double                                tally_multipl_factor = 1.0;
  bool                                  make_pwld = false;
  bool                                  uncollided_only = false;
  unsigned long long                    num_uncollided_particles = 1000;


  //derived from options-set during init
  bool                                  mesh_is_global = false;


//private:
//  chi_mesh::EmptyRegion empty_region;

public:
  //00
       Solver();
  //01
  bool Initialize();
  void InitMaterials();
  void InitParticleBatches();
  void InitTallies();
  void InitFieldFunctions();
  void InitGhostIDs();

  //02
  void Execute();
  void ExecuteRMCUncollided();

private:
  //03
  void Raytrace(Particle& prtcl);
  void RaytraceRMC(Particle& prtcl);
  void RaytraceUNC(Particle& prtcl);

  //04
  std::pair<int,chi_mesh::Vector3>
  ProcessScattering(Particle& prtcl,
                    chi_physics::TransportCrossSections* xs);
  void ProcessImportanceChange(Particle& prtcl);
  //05a
  void ContributeTally(Particle& prtcl,
                       chi_mesh::Vector3 pf);
  void ComputeRelativeStdDev();

  //05b
  double GetResidualFFPhi(std::vector<double>& N_in,
                          int dofs, int rmap,
                          chi_montecarlon::ResidualSource2* rsrc,
                          int egrp);
  chi_mesh::Vector3 GetResidualFFGradPhi(std::vector<chi_mesh::Vector3>& Grad_in,
                              int dofs, int rmap,
                              chi_montecarlon::ResidualSource2* rsrc,
                              int egrp);
  void ContributeTallyRMC(Particle& prtcl,
                          chi_mesh::Vector3 pf,
                          chi_mesh::RayDestinationInfo& ray_dest_info);

  void ContributeTallyUNC(Particle& prtcl,
                          chi_mesh::Vector3 pf,
                          double sig_t=0.0);

  //05c
  void RendesvouzTallies();
  void RendesvouzPWLTallies();

  //05d
  void NormalizeTallies();
  void NormalizePWLTallies();

  //06
  void ComputePWLDTransformations();
  void DevelopCollidedSource();

  //07
  void BuildMPITypes();
  void GetOutboundBankSize();
  void ReceiveIncomingParticles(std::vector<Particle>& inbound_particles);


  friend class chi_montecarlon::ResidualSource2;
  friend class chi_montecarlon::ResidualSource3;
};

#endif