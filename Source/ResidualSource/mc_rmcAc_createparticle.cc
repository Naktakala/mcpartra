#include "mc_rmcA_source.h"


#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>
#include "../../Solver/solver_montecarlon.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

////###################################################################
///**Executes a source sampling for the residual source.*/
//chi_montecarlon::Particle chi_montecarlon::ResidualSourceA::
//CreateParticle(chi_math::RandomNumberGenerator* rng)
//{
//  const double FOUR_PI = 4.0*M_PI;
//  chi_montecarlon::Particle new_particle;
//
//  std::vector<double>            shape_values;
//  std::vector<chi_mesh::Vector3> grad_shape_values;
//
//  //======================================== Choose interior or surface
//  bool sample_interior = false;
//  if (rng->Rand() < R_abs_localdomain_interior / (R_abs_localdomain_interior + R_abs_localdomain_surface))
//    sample_interior = true;
//
//  //################################################## INTERIOR
//  if (sample_interior)
//  {
//    int cell_local_id = std::lower_bound(
//      domain_cdf.begin(),
//      domain_cdf.end(),
//      rng->Rand()) - domain_cdf.begin();
//
//    //======================================== Randomly Sample Cell
//    auto& cell = ref_solver->grid->local_cells[cell_local_id];
//    auto cell_pwl_view =
//      ref_solver->pwl->MapFeViewL(cell_local_id);
//
//    //======================================== Sample energy group
//    int e_group = 0;
//    new_particle.egrp = e_group;
//
//    //======================================== Get material properties
//    int  mat_id         = cell.material_id;
//    int  xs_prop_id     = ref_solver->matid_xs_map[mat_id];
//    int  src_prop_id    = ref_solver->matid_q_map[mat_id];
//    auto material = chi_physics_handler.material_stack[mat_id];
//    auto xs = (chi_physics::TransportCrossSections*)material->properties[xs_prop_id];
//
//    double siga = xs->sigma_ag[e_group];
//    double Q    = 0.0;
//    if (src_prop_id >= 0)
//    {
//      auto prop = material->properties[src_prop_id];
//      auto q_prop = (chi_physics::IsotropicMultiGrpSource*)prop;
//      Q = q_prop->source_value_g[e_group];
//    }
//    //==================================== Sample position
//    chi_mesh::Vector3 pos =
//      GetRandomPositionInCell(*rng, cell_geometry_info[cell.local_id]);
//
//    new_particle.pos = pos;
//
//    //==================================== Sample direction
//    chi_mesh::Vector3 omega = RandomDirection(*rng);
//    new_particle.dir = omega;
//
//    //==================================== Populate shape values
//    cell_pwl_view->ShapeValues(pos, shape_values);
//    cell_pwl_view->GradShapeValues(pos,grad_shape_values);
//
//    //==================================== Get Residual
//    double phi = GetResidualFFPhi(shape_values,
//                                  cell_pwl_view->dofs,
//                                  cell.local_id,
//                                  e_group);
//    auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
//                                         cell_pwl_view->dofs,
//                                         cell.local_id,
//                                         e_group);
//
//    double r = (1.0/FOUR_PI)*
//               ( Q - siga*phi - omega.Dot(grad_phi) );
//
//    //======================================== Determine weight
//    if (std::fabs(r) > 1.0e-14)
//      new_particle.w = r / r_abs_cellk_interior_average[cell.local_id];
//    else
//    {
//      new_particle.w = 0.0;
//      new_particle.alive = false;
//    }
//
//    new_particle.cur_cell_local_id  = cell.local_id;
//    new_particle.cur_cell_global_id = cell.global_id;
//  }//interior
//  //################################################## SURFACE
//  else
//  {
//    int cf = std::lower_bound(
//      surface_cdf.begin(),
//      surface_cdf.end(),
//      rng->Rand()) - surface_cdf.begin();
//
//    auto& rcellface = r_abs_cellk_facef_surface_average[cf];
//
//    auto& cell = ref_solver->grid->local_cells[rcellface.cell_local_id];
//    auto cell_pwl_view =
//      ref_solver->pwl->MapFeViewL(rcellface.cell_local_id);
//
//    //==================================== Sample position
//    int f=rcellface.ass_face;
//    chi_mesh::Vector3 pos =
//      GetRandomPositionOnCellSurface(*rng,
//                                     cell_geometry_info[cell.local_id],
//                                     f);
//
//    auto& face = cell.faces[f];
//    auto& n = face.normal;
//
//    //======================================== Sample energy group
//    int e_group = 0;
//    new_particle.egrp = e_group;
//
//    new_particle.pos = pos;
//
//    //==================================== Sample direction
//    chi_mesh::Vector3 omega = RandomCosineLawDirection(*rng,-1.0*n);
//    new_particle.dir = omega;
//
//    //==================================== Populate shape values
//    cell_pwl_view->ShapeValues(pos, shape_values);
//
//    //==================================== Get Residual
//    double phi_P = GetResidualFFPhi(shape_values,
//                                    cell_pwl_view->dofs,
//                                    cell.local_id,
//                                    e_group);
//
//    double phi_N = phi_P;
//    if (not face.has_neighbor)
//      phi_N = 0.0; //TODO: Specialize for bndries
//
//    double r = (1.0/FOUR_PI)*(phi_N - phi_P);
//
//    //======================================== Determine weight
//    if (std::fabs(r) > 1.0e-14)
//      new_particle.w = r / rcellface.average_rstar;
//    else
//    {
//      new_particle.w = 0.0;
//      new_particle.alive = false;
//    }
//
//    new_particle.cur_cell_local_id  = cell.local_id;
//    new_particle.cur_cell_global_id = cell.global_id;
//
//  }//surface
//
//
//  return new_particle;
//}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSourceA::
CreateParticle(chi_math::RandomNumberGenerator* rng)
{
  const double FOUR_PI = 4.0*M_PI;
  chi_montecarlon::Particle new_particle;

  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  //======================================== Choose interior or surface
  bool sample_interior = false;
  if (rng->Rand() < R_abs_localdomain_interior / (R_abs_localdomain_interior + R_abs_localdomain_surface))
    sample_interior = true;

  //################################################## INTERIOR
  if (sample_interior)
  {
    //======================================== Randomly Sample Cell
    int cell_local_id = std::lower_bound(
      domain_cdf.begin(),
      domain_cdf.end(),
      rng->Rand()) - domain_cdf.begin();

    auto& cell = ref_solver->grid->local_cells[cell_local_id];
    auto& cell_pwl_view =
      ref_solver->pwl->GetCellFEView(cell_local_id);

    //======================================== Sample energy group
    int e_group = 0;
    new_particle.egrp = e_group;

    //======================================== Get material properties
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

    //======================================== Start rejection sampling
    bool particle_rejected = true;
    while (particle_rejected)
    {
      particle_rejected = false;
      new_particle.alive = true;

      //==================================== Sample position
      chi_mesh::Vector3 pos =
        GetRandomPositionInCell(*rng, cell_geometry_info[cell.local_id]);

      new_particle.pos = pos;

      //==================================== Sample direction
      chi_mesh::Vector3 omega = RandomDirection(*rng);
      new_particle.dir = omega;

      //==================================== Populate shape values
      cell_pwl_view.ShapeValues(pos, shape_values);
      cell_pwl_view.GradShapeValues(pos,grad_shape_values);

      //==================================== Get Residual
      double phi = GetResidualFFPhi(shape_values,
                                    cell_pwl_view.num_nodes,
                                    cell.local_id,
                                    e_group);
      auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
                                           cell_pwl_view.num_nodes,
                                           cell.local_id,
                                           e_group);

      double r = (1.0/FOUR_PI)*
                 ( Q - siga*phi - omega.Dot(grad_phi) );

      double rrandom = rng->Rand()*r_cellk_interior_max[cell.local_id];

      //======================================== Determine weight
      if (std::fabs(rrandom) < std::fabs(r))
        new_particle.w = r/std::fabs(r);
      else
      {
        new_particle.w = 0.0;
        new_particle.alive = false;
        particle_rejected = true;
      }

      new_particle.cur_cell_local_id  = cell.local_id;
      new_particle.cur_cell_global_id = cell.global_id;
    }//while particle rejected
  }//interior
    //################################################## SURFACE
  else
  {
    //==================================== Randomly sample face
    int cf = std::lower_bound(
      surface_cdf.begin(),
      surface_cdf.end(),
      rng->Rand()) - surface_cdf.begin();

    auto& rcellface = r_abs_cellk_facef_surface_average[cf];

    auto& cell = ref_solver->grid->local_cells[rcellface.cell_local_id];
    auto& cell_pwl_view =
      ref_solver->pwl->GetCellFEView(rcellface.cell_local_id);

    //======================================== Start rejection sampling
    bool particle_rejected = true;
    while (particle_rejected)
    {
      particle_rejected = false;
      new_particle.alive = true;

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
      cell_pwl_view.ShapeValues(pos, shape_values);

      //==================================== Get Residual
      double phi_P = GetResidualFFPhi(shape_values,
                                      cell_pwl_view.num_nodes,
                                      cell.local_id,
                                      e_group);

      double phi_N = phi_P;
      if (not face.has_neighbor)
        phi_N = 0.0; //TODO: Specialize for bndries

      double r = (1.0/FOUR_PI)*(phi_N - phi_P);

      double rrandom = rng->Rand()*rcellface.maximum;

      //======================================== Determine weight
      if (std::fabs(rrandom) < std::fabs(r))
        new_particle.w = r/std::fabs(r);
      else
      {
        new_particle.w = 0.0;
        new_particle.alive = false;
        particle_rejected = true;
      }

      new_particle.cur_cell_local_id  = cell.local_id;
      new_particle.cur_cell_global_id = cell.global_id;
    }//while particle rejected
  }//surface


  return new_particle;
}

double chi_montecarlon::ResidualSourceA::
  GetParallelRelativeSourceWeight()
{
  double global_total_source_weight =
    (R_abs_globaldomain_interior + R_abs_globaldomain_surface);

  relative_weight = (R_abs_localdomain_interior + R_abs_localdomain_surface) / global_total_source_weight;
  return  relative_weight;
}