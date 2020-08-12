#include "mc_rmc3_source.h"

#include "../../RandomNumberGenerator/montecarlon_rng.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>
#include "../../Solver/solver_montecarlon.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualSource3::
ResidualSource3(chi_physics::FieldFunction *in_resid_ff,
                bool use_uniform_sampling) :
  sample_uniformly(use_uniform_sampling)
{
  type_index = SourceTypes::RESIDUAL3;
  resid_ff = in_resid_ff;
  particles_L = 0;
  particles_R = 0;

  weights_L = 0.0;
  weights_R = 0.0;
}

//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void chi_montecarlon::ResidualSource3::
Initialize(chi_mesh::MeshContinuum *ref_grid,
           SpatialDiscretization_FV *ref_fv_sdm,
           chi_montecarlon::Solver* ref_solver)
{
  const double FOUR_PI = 4.0*M_PI;
  chi_log.Log(LOG_0) << "Initializing Residual3 Sources";
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  this->ref_solver = ref_solver;
  RandomNumberGenerator& rng = ref_solver->rng0;
  size_t num_local_cells = ref_solver->grid->local_cells.size();

  BuildCellVolInfo(ref_grid,ref_fv_sdm);

  RemoveFFDiscontinuities();

  //============================================= Sample each cell
  chi_log.Log(LOG_0) << "Integrating cell source.";
  cell_avg_interior_rstar.clear();
  cell_avg_interior_rstar.resize(num_local_cells, 0.0);
  cell_avg_surface_rstar.clear();
  cell_avg_surface_rstar.resize(num_local_cells, 0.0);
  cell_volumes.clear();
  cell_volumes.resize(num_local_cells,0.0);
  cell_IntVOmega_rstar.clear();
  cell_IntVOmega_rstar.resize(num_local_cells, 0.0);
  cell_IntSOmega_rstar.clear();
  cell_IntSOmega_rstar.resize(num_local_cells, 0.0);

  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;
  for (auto& cell : ref_solver->grid->local_cells)
  {
    chi_log.Log(LOG_0) << "Cell " << cell.local_id;
    auto cell_pwl_view  =
      ref_solver->pwl_discretization->MapFeViewL(cell.local_id);
    auto cell_FV_view =
      ref_solver->fv_discretization->MapFeView(cell.local_id);

    int e_group = 0;

    //====================================== Get material properties
    int  mat_id         = cell.material_id;
    int  xs_prop_id     = ref_solver->matid_xs_map[mat_id];
    int  src_prop_id    = ref_solver->matid_q_map[mat_id];
    auto material = chi_physics_handler.material_stack[mat_id];
    auto xs = (chi_physics::TransportCrossSections*)material->properties[xs_prop_id];

    double siga = xs->sigma_ag[e_group];
    double Q    = 0.0;
    if (src_prop_id >= 0)
    {
      auto prop = material->properties[src_prop_id];
      auto q_prop = (chi_physics::IsotropicMultiGrpSource*)prop;
      Q = q_prop->source_value_g[e_group];
    }

    //====================================== Sample points
    int num_points = 1000;
    for (int k=0; k<num_points; ++k)
    {
      //==================================== Sample direction
      chi_mesh::Vector3 omega = RandomDirection(rng);

      //==================================== Sample position
      chi_mesh::Vector3 pos =
        GetRandomPositionInCell(rng,
                                cell_geometry_info[cell.local_id]);

      //==================================== Populate shape values
      cell_pwl_view->ShapeValues(pos, shape_values);
      cell_pwl_view->GradShapeValues(pos,grad_shape_values);

      //==================================== Get Residual
      double phi = GetResidualFFPhi(shape_values,
                                    cell_pwl_view->dofs,
                                    cell.local_id,
                                    e_group);
      auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
                                           cell_pwl_view->dofs,
                                           cell.local_id,
                                           e_group);
      double r = (1.0/FOUR_PI)*
                 ( Q - siga*phi - omega.Dot(grad_phi) );

      //==================================== Contribute to avg
      cell_avg_interior_rstar[cell.local_id] += std::fabs(r);
    }
    cell_avg_interior_rstar[cell.local_id]  /= num_points;
    cell_volumes[cell.local_id]    = cell_FV_view->volume;
    cell_IntVOmega_rstar[cell.local_id] =
      cell_volumes[cell.local_id] * FOUR_PI *
      cell_avg_interior_rstar[cell.local_id];

    domain_volume += cell_FV_view->volume;
    if (std::fabs(cell_avg_interior_rstar[cell.local_id]) > 1.0e-16)
      source_volume += cell_FV_view->volume;
  }//for cell

  MPI_Barrier(MPI_COMM_WORLD);

  //============================================= Sample surfaces
  chi_log.Log(LOG_0) << "Integrating surface source.";
  for (auto& cell : ref_solver->grid->local_cells)
  {
    auto cell_pwl_view =
      ref_solver->pwl_discretization->MapFeViewL(cell.local_id);
    auto cell_FV_view =
      ref_solver->fv_discretization->MapFeView(cell.local_id);

    int e_group = 0;

    int f=-1;
    for (auto& face : cell.faces)
    {
      ++f;
      double A_f = cell_FV_view->face_area[f];
      auto&  n   = face.normal;

      RCellFace rcellface;
      rcellface.cell_local_id = cell.local_id;
      rcellface.ass_face      = f;
      rcellface.area          = A_f;

      double face_ave_surface_pstar=0.0;
      int num_points = 250;
      for (int k=0; k<num_points; ++k)
      {
        //==================================== Sample position
        chi_mesh::Vector3 pos =
          GetRandomPositionOnCellSurface(rng,
                                         cell_geometry_info[cell.local_id],
                                         f);

        //==================================== Populate shape values
        cell_pwl_view->ShapeValues(pos, shape_values);

        //==================================== Get Residual
        double phi_P = GetResidualFFPhi(shape_values,
                                        cell_pwl_view->dofs,
                                        cell.local_id,
                                        e_group);

        double phi_N = phi_P;
        if (face.neighbor < 0)
          phi_N = 0.0; //TODO: Specialize for bndries
        else
        {
//          int adj_local_id = face.GetNeighborLocalID(grid);
//          auto neighbor_pwl_view =
//            ref_solver->pwl_discretization->MapFeViewL(adj_local_id);
//
//          neighbor_pwl_view->ShapeValues(pos, shape_values);
//          phi_N = GetResidualFFPhi(shape_values,
//                                   neighbor_pwl_view->dofs,
//                                   adj_local_id,
//                                   e_group);
        }

        double r = (1.0/FOUR_PI)*(phi_N - phi_P);

        //==================================== Contribute to avg
        face_ave_surface_pstar += std::fabs(r);
      }//for k
      face_ave_surface_pstar /= num_points;

      rcellface.average_rstar = face_ave_surface_pstar;

      cell_face_pairs.push_back(rcellface);
    }//for face
  }//for cell

  //============================================= Integrate sources
  IntVOmega_rstar = 0.0;
  for (double v : cell_IntVOmega_rstar)
    IntVOmega_rstar += v;

  chi_log.Log(LOG_ALL) << "Total interior source: " << IntVOmega_rstar;

  IntSOmega_rstar = 0.0;
  for (auto& rcellface : cell_face_pairs)
    IntSOmega_rstar += rcellface.average_rstar*rcellface.area*M_PI;

  chi_log.Log(LOG_ALL) << "Total surface source: " << IntSOmega_rstar;

  //============================================= Compute interior cdf
  {
    interior_cdf.clear();
    interior_cdf.resize(num_local_cells, 0.0);
    double running_total = 0.0;
    for (int c = 0; c < num_local_cells; ++c)
    {
      running_total += cell_IntVOmega_rstar[c];
      interior_cdf[c] = running_total / IntVOmega_rstar;
    }
  }

  //============================================= Compute surface cdf
  {
    surface_cdf.clear();
    surface_cdf.resize(cell_face_pairs.size(), 0.0);
    double running_total = 0.0;
    for (int cf = 0; cf < cell_face_pairs.size(); ++cf)
    {
      auto& rcellface = cell_face_pairs[cf];
      running_total += rcellface.average_rstar*rcellface.area*M_PI;
      surface_cdf[cf] = running_total / IntSOmega_rstar;
    }
  }
  chi_log.Log(LOG_0) << "Done initializing Residual Sources";

  MPI_Allreduce(&IntVOmega_rstar,
                &IntVOmega_rstar_global,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  MPI_Allreduce(&IntSOmega_rstar,
                &IntSOmega_rstar_global,
                1,
                MPI_DOUBLE,
                MPI_SUM,
                MPI_COMM_WORLD);

  chi_log.Log(LOG_0) << "Total interior source: " << IntVOmega_rstar_global;
  chi_log.Log(LOG_0) << "Total surface source: " << IntSOmega_rstar_global;

}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource3::
CreateParticle(chi_montecarlon::RandomNumberGenerator* rng)
{
  const double FOUR_PI = 4.0*M_PI;
  chi_montecarlon::Particle new_particle;

  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  //======================================== Choose interior or surface
  bool sample_interior = false;
  if (rng->Rand() < IntVOmega_rstar / (IntVOmega_rstar + IntSOmega_rstar))
    sample_interior = true;

  //################################################## INTERIOR
  if (sample_interior)
  {
    int cell_local_id = std::lower_bound(
      interior_cdf.begin(),
      interior_cdf.end(),
      rng->Rand()) - interior_cdf.begin();

    //======================================== Randomly Sample Cell
    auto& cell = ref_solver->grid->local_cells[cell_local_id];
    auto cell_pwl_view =
      ref_solver->pwl_discretization->MapFeViewL(cell_local_id);

    //======================================== Sample energy group
    int e_group = 0;
    new_particle.egrp = e_group;

    //======================================== Get material properties
    int  mat_id         = cell.material_id;
    int  xs_prop_id     = ref_solver->matid_xs_map[mat_id];
    int  src_prop_id    = ref_solver->matid_q_map[mat_id];
    auto material = chi_physics_handler.material_stack[mat_id];
    auto xs = (chi_physics::TransportCrossSections*)material->properties[xs_prop_id];

    double siga = xs->sigma_tg[e_group];
    double Q    = 0.0;
    if (src_prop_id >= 0)
    {
      auto prop = material->properties[src_prop_id];
      auto q_prop = (chi_physics::IsotropicMultiGrpSource*)prop;
      Q = q_prop->source_value_g[e_group];
    }
    //==================================== Sample position
    chi_mesh::Vector3 pos =
      GetRandomPositionInCell(*rng, cell_geometry_info[cell.local_id]);

    new_particle.pos = pos;

    //==================================== Sample direction
    chi_mesh::Vector3 omega = RandomDirection(*rng);
    new_particle.dir = omega;

    //==================================== Populate shape values
    cell_pwl_view->ShapeValues(pos, shape_values);
    cell_pwl_view->GradShapeValues(pos,grad_shape_values);

    //==================================== Get Residual
    double phi = GetResidualFFPhi(shape_values,
                                  cell_pwl_view->dofs,
                                  cell.local_id,
                                  e_group);
    auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
                                         cell_pwl_view->dofs,
                                         cell.local_id,
                                         e_group);

    double r = (1.0/FOUR_PI)*
               ( Q - siga*phi - omega.Dot(grad_phi) );

    //======================================== Determine weight
    if (std::fabs(r) > 1.0e-14)
      new_particle.w = r * (IntVOmega_rstar_global + IntSOmega_rstar_global) /
                       cell_avg_interior_rstar[cell.local_id];
    else
    {
      new_particle.w = 0.0;
      new_particle.alive = false;
    }
//  new_particle.alive = false;

    new_particle.cur_cell_local_id  = cell.local_id;
    new_particle.cur_cell_global_id = cell.global_id;
  }//interior
  //################################################## SURFACE
  else
  {
    int cf = std::lower_bound(
      surface_cdf.begin(),
      surface_cdf.end(),
      rng->Rand()) - surface_cdf.begin();

    auto& rcellface = cell_face_pairs[cf];

    auto& cell = ref_solver->grid->local_cells[rcellface.cell_local_id];
    auto cell_pwl_view =
      ref_solver->pwl_discretization->MapFeViewL(rcellface.cell_local_id);

    //================================================ Sample face
    //==================================== Sample position
    int f=rcellface.ass_face;
    chi_mesh::Vector3 pos =
      GetRandomPositionOnCellSurface(*rng,
                                     cell_geometry_info[cell.local_id],
                                     f);

    auto& face = cell.faces[f];
    auto& n = face.normal;

    //======================================== Sample energy group
    int e_group = 0;
    new_particle.egrp = e_group;

    new_particle.pos = pos;

    //==================================== Sample direction
    chi_mesh::Vector3 omega = RandomCosineLawDirection(*rng,-1.0*n);
    new_particle.dir = omega;

    //==================================== Populate shape values
    cell_pwl_view->ShapeValues(pos, shape_values);

    //==================================== Get Residual
    double phi_P = GetResidualFFPhi(shape_values,
                                    cell_pwl_view->dofs,
                                    cell.local_id,
                                    e_group);

    double phi_N = phi_P;
    if (face.neighbor < 0)
      phi_N = 0.0; //TODO: Specialize for bndries
    else
    {
//      int adj_local_id = face.GetNeighborLocalID(grid);
//      auto neighbor_pwl_view =
//        ref_solver->pwl_discretization->MapFeViewL(adj_local_id);
//
//      neighbor_pwl_view->ShapeValues(pos, shape_values);
//      phi_N = GetResidualFFPhi(shape_values,
//                               neighbor_pwl_view->dofs,
//                               adj_local_id,
//                               e_group);
    }

    double r = (1.0/FOUR_PI)*(phi_N - phi_P);

    //======================================== Determine weight
    if (std::fabs(r) > 1.0e-14)
      new_particle.w = r * (IntVOmega_rstar_global + IntSOmega_rstar_global) /
                       rcellface.average_rstar;
    else
    {
      new_particle.w = 0.0;
      new_particle.alive = false;
    }
    //  new_particle.alive = false;

    new_particle.cur_cell_local_id  = cell.local_id;
    new_particle.cur_cell_global_id = cell.global_id;

  }//surface


  return new_particle;
}

double chi_montecarlon::ResidualSource3::
  GetParallelRelativeSourceWeight()
{
  double global_total_source_weight =
    (IntVOmega_rstar_global + IntSOmega_rstar_global);

  relative_weight = (IntVOmega_rstar + IntSOmega_rstar)/global_total_source_weight;
  return  relative_weight;
}