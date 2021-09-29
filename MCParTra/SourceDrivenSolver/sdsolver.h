#ifndef SSSOLVER_MONTECARLON_H
#define SSSOLVER_MONTECARLON_H

#include "mcpartra.h"
#include "mcpartra_particle.h"
#include "Utils/MultigroupTally/multigroup_tally.h"
#include "Utils/CustomVolumeTally/customvolumetally.h"
#include "Raytracing/raytracing.h"
#include "Sources/mc_base_source.h"
#include "Sources/mc_volume_src_element.h"
#include "Sources/mc_surface_src_element.h"
#include "MGScatteringCDFs/mg_scattering_cdfs.h"

#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"
#include "ChiMesh/Raytrace/raytracing.h"

#include "ChiMath/chi_math.h"
#include "ChiMath/RandomNumberGeneration/random_number_generator.h"
#include "ChiMath/UnknownManager/unknown_manager.h"

#include "ChiPhysics/SolverBase/chi_solver.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"
#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include "ChiMath/SpatialDiscretization/FiniteVolume/fv.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>

namespace mcpartra
{

//######################################################### Class def
/**Monte Carlo neutron particle solver.*/
class SourceDrivenSolver : public chi_physics::Solver
{
protected:
  typedef std::shared_ptr<SpatialDiscretization_FV> SDMFVPtr;
  typedef std::shared_ptr<SpatialDiscretization_PWLD> SDMPWLPtr;
public:

  /**Bit-wise identifier to a specific tally. This allows
   * particles to ship with a bit-wise integer that functions
   * as 8bits*4bytes=32 flag values.*/
  enum TallyMask : unsigned int
  {
    DEFAULT_FVTALLY       = 1 << 0, //0000 0001
    DEFAULT_PWLTALLY      = 1 << 1, //0000 0010
    UNCOLLIDED_FVTALLY    = 1 << 2, //0000 0100
    UNCOLLIDED_PWLTALLY   = 1 << 3, //0000 1000
  };

  /**Maps a bit-wise tally identifier to an actual index.*/
  std::map<TallyMask,unsigned int> TallyMaskIndex =
    {{DEFAULT_FVTALLY,     0},
     {DEFAULT_PWLTALLY,    1},
     {UNCOLLIDED_FVTALLY,  2},
     {UNCOLLIDED_PWLTALLY, 3}};

  std::vector<unsigned int> fv_tallies  = {TallyMaskIndex[DEFAULT_FVTALLY]};
  std::vector<unsigned int> pwl_tallies = {TallyMaskIndex[DEFAULT_PWLTALLY]};

  //=================================== Members
private:
  chi_mesh::MeshContinuumPtr            grid = nullptr;

  // Materials related members
  std::map<int, std::shared_ptr<chi_physics::TransportCrossSections>> matid_xs_map2;
  std::vector<bool>                       matid_has_q_flags;
  std::map<int, std::shared_ptr<chi_physics::IsotropicMultiGrpSource>> matid_q_map2;
  std::map<int, MultigroupScatteringCDFs> matid_scattering_cdfs;

  // Tally related members
  SDMFVPtr                              fv;
  SDMPWLPtr                             pwl;
  chi_math::UnknownManager              uk_man_fv;
  chi_math::UnknownManager              uk_man_pwld;
  std::vector<MultigroupTally>          grid_tally_blocks;
  std::vector<CustomVolumeTally>        custom_tallies;

  // Variance reduction quantities
public:
  chi_math::UnknownManager              uk_man_fv_importance;
  size_t                                importance_num_groups=0;
  std::vector<double>                   local_cell_importance_setting;
  std::vector<chi_mesh::Vector3>        local_cell_importance_directions;
  std::vector<std::pair<double,double>> local_cell_importance_exp_coeffs;
public:
  std::vector<double>                   local_cell_importance;

  // Source information
public:
  std::vector<CellGeometryData>         cell_geometry_info;
  std::vector<SourceBase*>              sources;

private:
  double                                accumulated_src_importances = 0.0;
  double                                total_local_source_rate = 0.0;
  double                                total_globl_source_rate = 0.0;
  std::vector<double>                   local_source_cdf;
  double                                source_normalization = 1.0;

  // Particle batch information
  std::vector<unsigned long long>       batch_sizes;
  std::vector<unsigned long long>       batch_sizes_per_loc;

  std::vector<unsigned long long>       uncollided_batch_sizes;
  std::vector<unsigned long long>       uncollided_batch_sizes_per_loc;

public:
  size_t                                num_groups=1; //updated during material init
  size_t                                num_moments=1;

public:
  chi_math::RandomNumberGenerator       rng0;
  chi_mesh::RayTracer                   default_raytracer;


  //========================= Runtime quantities
private:
  std::vector<double>                   segment_lengths;
  std::vector<double>                   N_f, N_i;
  unsigned long long                    nps=0;
  unsigned long long                    nps_global=0;

  double                                max_sigma=0.0;
  double                                max_relative_sigma=0.0;
  double                                avg_sigma=0.0;

  double                                max_fem_sigma=0.0;

  MPI_Datatype                          mpi_prtcl_data_type=0;
  std::vector<Particle>                 outbound_particle_bank;
  std::vector<Particle>                 particle_source_bank;
  unsigned int                          total_outbound_bank_size=0;

  std::vector<std::string>              lost_particles;
public:
  std::vector<std::pair<int,int>>       m_to_ell_em_map;

public:


  //========================= Options
  struct Options
  {
    unsigned long long num_particles = 1000;
    unsigned long long num_uncollided_particles = 1000;
    unsigned long long tally_rendezvous_intvl = 100000;
    bool               mono_energy = false;
    int                scattering_order = 1; ///< Maximum 7
    bool               force_isotropic = false;
    int                group_hi_bound = -1;
    int                group_lo_bound = -1;
    double             tally_multipl_factor = 1.0;
    bool               make_pwld = false;
    bool               uncollided_only = false;
    bool               importances_during_raytracing = false;
    bool               apply_source_importance_sampling = false;

    bool               write_run_tape = false;
    std::string        run_tape_base_name;
  }options;

public:
  //00
  explicit SourceDrivenSolver(const std::string& text_name);
  explicit SourceDrivenSolver(int seed, const std::string& text_name);

  //01
  void Initialize() override;
  void InitRaytracing();
  void InitMaterials(); //01a
  void InitCellImportances(); //01b
  void InitMomentIndices(); //01d
  void InitTallies(); //01e
  void InitFieldFunctions(); //01f
  void InitCellGeometryData();
  void InitSources(); //01h
  void InitParticleBatches(); //01i

  //02
  void Execute() override;
  void DepleteSourceParticleBank();
  void CommAndSimOutboundInboundParticles();
  //02a
  Particle SampleSources(chi_math::RandomNumberGenerator& rng);
  //02b
  void PrintBatchInfo(size_t b, double particle_rate);

private:
  //03
  void Raytrace(Particle& prtcl);

  //04
  std::pair<int,chi_mesh::Vector3>
  ProcessScattering(Particle& prtcl,
                    std::shared_ptr<chi_physics::TransportCrossSections> xs);
  std::pair<int,chi_mesh::Vector3>
  ProcessScattering(Particle& prtcl,
                    const MultigroupScatteringCDFs& xs);
  static std::pair<double, double>
  GetMaterialSigmas(const chi_physics::TransportCrossSections& xs,
                    unsigned int energy_group) ;

  void ProcessImportanceChange(Particle& prtcl, double current_cell_importance);

  //05a
  void ContributeTally(const chi_mesh::Cell& cell,
                       const Particle& prtcl,
                       const chi_mesh::Vector3& pf);

  double ContributeTallyUNC(const chi_mesh::Cell& cell,
                            const Particle& prtcl,
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
//  void DevelopCollidedSource(chi_montecarlon::MultigroupTally& input_fv_tally);

  //General utils
public:
  size_t AddCustomVolumeTally(chi_mesh::LogicalVolume& logical_volume);

  //IO Utils
  void WriteRunTape(const std::string& file_base_name);
  void ReadRunTape(const std::string& file_name);
  void WriteLBSFluxMoments(const std::string& file_name);
  void ReadImportanceMap(const std::string& file_name);

  friend class ResidualSourceB;
  friend class ResidualSourceA;

};
}//namespace mcpartra
#endif //SSSOLVER_MONTECARLON_H