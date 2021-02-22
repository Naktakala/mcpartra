#ifndef SOLVER_MONTECARLON_H
#define SOLVER_MONTECARLON_H

#include"../chi_montecarlon.h"
#include "ChiMath/RandomNumberGeneration/random_number_generator.h"
#include "../chi_montecarlon_particle.h"
#include"ChiPhysics/SolverBase/chi_solver.h"
#include "../Source/mc_base_source.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiPhysics/PhysicsMaterial/property10_transportxsections.h"
#include "FiniteVolume/fv.h"
#include "FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/chi_math.h"

#include "ChiMesh/Raytrace/raytracing.h"

#include "ChiMath/UnknownManager/unknown_manager.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>

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

  class GridTallyBlock;
  class CustomVolumeTally;
}

//######################################################### Tally struct
/**Tally block.*/
class chi_montecarlon::GridTallyBlock
{
private:
  bool is_empty = true;
public:
  std::vector<double> tally_local;
  std::vector<double> tally_global;
  std::vector<double> tally_sqr_local;
  std::vector<double> tally_sqr_global;

  std::vector<double> tally_sigma;
  std::vector<double> tally_relative_sigma;

  void Resize(int tally_size)
  {
    tally_local         .resize(tally_size, 0.0);
    tally_global        .resize(tally_size, 0.0);
    tally_sqr_local     .resize(tally_size, 0.0);
    tally_sqr_global    .resize(tally_size, 0.0);
    tally_sigma         .resize(tally_size, 0.0);
    tally_relative_sigma.resize(tally_size, 0.0);

    is_empty = false;
  }

  GridTallyBlock& operator=(const GridTallyBlock& that)
  {
    tally_local          = that.tally_local;
    tally_global         = that.tally_global;
    tally_sqr_local      = that.tally_sqr_local;
    tally_sqr_global     = that.tally_sqr_global;

    tally_sigma          = that.tally_sigma;
    tally_relative_sigma = that.tally_relative_sigma;

    return *this;
  }

  void ZeroOut()
  {
    int tally_size = tally_local.size();

    tally_local          .assign(tally_size, 0.0);
    tally_global         .assign(tally_size, 0.0);
    tally_sqr_local      .assign(tally_size, 0.0);
    tally_sqr_global     .assign(tally_size, 0.0);

    tally_sigma          .assign(tally_size, 0.0);
    tally_relative_sigma .assign(tally_size, 0.0);
  }

  GridTallyBlock& operator+=(const GridTallyBlock& that)
  {
    int tally_size = tally_local.size();

    for (int i=0; i<tally_size; ++i)
    {
      tally_local         [i] += that.tally_local         [i];
      tally_global        [i] += that.tally_global        [i];
      tally_sqr_local     [i] += that.tally_sqr_local     [i];
      tally_sqr_global    [i] += that.tally_sqr_global    [i];
      tally_sigma         [i] += that.tally_sigma         [i];
      tally_relative_sigma[i] += that.tally_relative_sigma[i];
    }

    return *this;
  }

  GridTallyBlock& operator*=(const double value)
  {
    int tally_size = tally_local.size();

    for (int i=0; i<tally_size; ++i)
    {
      tally_local         [i] *= value;
      tally_global        [i] *= value;
      tally_sqr_local     [i] *= value;
      tally_sqr_global    [i] *= value;
      tally_sigma         [i] *= value;
      tally_relative_sigma[i] *= value;
    }

    return *this;
  }

  bool empty() {return is_empty;}
};

//#########################################################
/**Custom tally Structure.*/
class chi_montecarlon::CustomVolumeTally
{
public:
  struct TallyFluctuationChart
  {
    std::vector<double> average;
    std::vector<double> sigma;
  };
  typedef TallyFluctuationChart TFC;
public:
  GridTallyBlock        grid_tally;
  std::vector<bool>     local_cell_tally_mask; ///< Indicates whether a local cell is part of tally
  double                tally_volume=0.0;
  bool                  initialized=false;
  std::vector<TFC>      tally_fluctation_chart;

  explicit
  CustomVolumeTally(std::vector<bool>& in_masking) :
    local_cell_tally_mask(in_masking) {}

  void Initialize(size_t tally_size, double in_volume)
  {
    grid_tally.Resize(tally_size);
    tally_volume = in_volume;
    initialized = true;
  }
};

//######################################################### Class def
/**Monte Carlo neutron particle solver.*/
class chi_montecarlon::Solver : public chi_physics::Solver
{
protected:
  typedef std::shared_ptr<SpatialDiscretization_FV> SDMFVPtr;
  typedef std::shared_ptr<SpatialDiscretization_PWL> SDMPWLPtr;
public:
  enum RayTraceMethod
  {
    STANDARD   = 0,
    UNCOLLIDED = 1
  };

  struct TallyMethod
  {
    static const int STANDARD     = 0;
    static const int RMC_CHAR_RAY = 1;
  };

  /**Bit-wise identifier to a specific tally. This allows
   * particles to ship with a bit-wise integer that functions
   * as 8bits*4bytes=32 flag values.*/
  enum TallyMask
  {
    DEFAULT_FVTALLY       = 1 << 0, //0000 0001
    DEFAULT_PWLTALLY      = 1 << 1, //0000 0010
    UNCOLLIDED_FVTALLY    = 1 << 2, //0000 0100
    UNCOLLIDED_PWLTALLY   = 1 << 3, //0000 1000
    MAKE_DIRECT_PARTICLES = 1 << 4, //0001 0000
  };

  /**Maps a bit-wise tally identifier to an actual index.*/
  std::map<TallyMask,int> TallyMaskIndex =
    {{DEFAULT_FVTALLY,     0},
     {DEFAULT_PWLTALLY,    1},
     {UNCOLLIDED_FVTALLY,  2},
     {UNCOLLIDED_PWLTALLY, 3}};

  std::vector<int> fv_tallies =
    {TallyMaskIndex[DEFAULT_FVTALLY]};

  std::vector<int> pwl_tallies =
    {TallyMaskIndex[DEFAULT_PWLTALLY]};

private:
  chi_mesh::MeshContinuumPtr            grid;
  std::map<int,int>                     cell_neighbor_nonlocal_local_id;
  SDMFVPtr                              fv;
  SDMPWLPtr                             pwl;
private:
  std::vector<int>                      matid_xs_map;
  std::vector<int>                      matid_q_map;

  std::vector<unsigned long long>       batch_sizes;
  std::vector<unsigned long long>       batch_sizes_per_loc;

  std::vector<unsigned long long>       uncollided_batch_sizes;
  std::vector<unsigned long long>       uncollided_batch_sizes_per_loc;
public:
  int                                   num_grps=1; //updated during material init
  int                                   num_moms=1;
  chi_math::UnknownManager              uk_man_fv;
  chi_math::UnknownManager              uk_man_pwld;

public:
  std::vector<GridTallyBlock>           grid_tally_blocks;
  std::vector<CustomVolumeTally>        custom_tallies;
  double                                source_normalization = 1.0;

public:
  chi_math::RandomNumberGenerator       rng0;
  std::vector<chi_montecarlon::Source*> sources;

  //======================== Variance reduction quantities
public:
  std::vector<double>                   local_cell_importance_setting;
private:
  std::vector<double>                   local_cell_importance;

  //======================== RMC quantities
public:
  std::vector<double>                   cdf_phi_unc_group;
  std::vector<std::vector<double>>      cdf_phi_unc_group_cell;
  std::vector<std::vector<double>>      IntVk_phi_unc_g;
  double                                IntVSumG_phi_unc=0.0;
  std::vector<double>                   IntV_phi_unc_g;
  double                                domain_volume = 0.0;
  //========================= Runtime quantities
private:
  std::vector<double>                   segment_lengths;
  std::vector<double>                   N_f,N_i;
  std::vector<chi_mesh::Vector3>        Grad;
  size_t                                current_batch;
  unsigned long long                    nps;
  unsigned long long                    nps_global;

  double                                max_sigma=0.0;
  double                                max_relative_sigma=0.0;
  double                                avg_sigma=0.0;

  double                                max_fem_sigma=0.0;
  double                                max_fem_relative_sigma=0.0;
  double                                avg_fem_sigma=0.0;

  MPI_Datatype                          mpi_prtcl_data_type;
  std::vector<Particle>                 outbound_particle_bank;
  std::vector<Particle>                 particle_source_bank;
  unsigned int                          total_outbound_bank_size=0;

public:


  //========================= Options
  double                                tolerance = 0.000001;
  unsigned long long                    num_particles = 1000;
  unsigned long long                    num_uncollided_particles = 1000;
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


  //derived from options-set during init
  bool                                  mesh_is_global = false;

public:
  //00
       Solver();
       std::pair<int,int> GetGroupBounds();
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
  void PrintBatchInfo(int b, double particle_rate);

private:
  //03
  void Raytrace(Particle& prtcl);
  void RaytraceSTD(Particle& prtcl);
  void RaytraceUNC(Particle& prtcl);
  void RaytraceRMC(Particle& prtcl);

  //04
  std::pair<int,chi_mesh::Vector3>
  ProcessScattering(Particle& prtcl,
                    chi_physics::TransportCrossSections* xs);
  void ProcessImportanceChange(Particle& prtcl);

  //05a
  void ContributeTally(Particle& prtcl,
                       const chi_mesh::Vector3& pf);
  double
  GetResidualFFPhi(std::vector<double>& N_in,
                   int dofs, int rmap,
                   chi_montecarlon::ResidualSourceB* rsrc,
                   int egrp);
  chi_mesh::Vector3
  GetResidualFFGradPhi(std::vector<chi_mesh::Vector3>& Grad_in,
                       int dofs, int rmap,
                       chi_montecarlon::ResidualSourceB* rsrc,
                       int egrp);
  void ContributeTallyRMC(Particle& prtcl,
                          const chi_mesh::Vector3& pf,
                          chi_mesh::RayDestinationInfo& ray_dest_info);
  Particle MakeScatteredParticle(Particle& prtcl,
                             double tracklength,
                             double weight);

  void ContributeTallyUNC(Particle& prtcl,
                          const chi_mesh::Vector3& pf,
                          double sig_t=0.0);

  //05b
  void RendesvouzTallies();

  //05c
  void ComputeUncertainty();

  //05d
  void NormalizeTallies();

  //06
  void ComputePWLDTransformations();

  //07
  void BuildMPITypes();
  void GetOutboundBankSize();
  void ReceiveIncomingParticles(std::vector<Particle>& inbound_particles);

  //08
  void DevelopCollidedSource(chi_montecarlon::GridTallyBlock& input_fv_tally);

  //General utils
public:
  size_t AddCustomVolumeTally(chi_mesh::LogicalVolume& logical_volume);

  friend class chi_montecarlon::ResidualSourceB;
  friend class chi_montecarlon::ResidualSourceA;

};

#endif