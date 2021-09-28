#include "mc_rmcA_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Executes a source sampling for the residual source.*/
mcpartra::Particle mcpartra::ResidualSourceA::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  constexpr double FOUR_PI = 4.0*M_PI;
  mcpartra::Particle new_particle;

  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  //======================================== Sample energy group
  int e_group = 0;
  new_particle.egrp = e_group;

  //======================================== Choose interior or surface
  bool sample_interior = false;
  if (rng.Rand() < R_abs_localdomain_interior /
          (R_abs_localdomain_interior + R_abs_localdomain_surface))
  {sample_interior = true;}

  //################################################## INTERIOR
  if (sample_interior)
  {
    //======================================== Randomly Sample Cell
    int64_t cell_local_id = std::lower_bound(domain_cdf.begin(),
                                             domain_cdf.end(),
                                             rng.Rand()) - domain_cdf.begin();

    const auto&  cell           = ref_solver.grid->local_cells[cell_local_id];
    const auto&  cell_pwl_view  = ref_solver.pwl->GetCellMappingFE(cell_local_id);
    const auto&  cell_geom_info = cell_geometry_info->operator[](cell_local_id);
    const size_t cell_num_nodes = ref_solver.pwl->GetCellNumNodes(cell);
    const auto&  cell_r_info    = residual_info_cell_interiors[cell_local_id];

    //====================================== Get material properties
    MaterialData mat_data;
    PopulateMaterialData(cell.material_id,e_group,mat_data);

    auto siga = mat_data.siga;
    auto Q    = mat_data.Q;

    //======================================== Start rejection sampling
    bool particle_rejected = true;
    while (particle_rejected)
    {
      particle_rejected = false;
      new_particle.alive = true;

      //==================================== Sample position
      new_particle.pos = GetRandomPositionInCell(rng, cell_geom_info);

      //==================================== Sample direction
      chi_mesh::Vector3 omega = SampleRandomDirection(rng);
      new_particle.dir = omega;

      //==================================== Populate shape values
      cell_pwl_view->ShapeValues(new_particle.pos, shape_values);
      cell_pwl_view->GradShapeValues(new_particle.pos, grad_shape_values);

      //==================================== Get Residual
      double phi = GetResidualFFPhi(shape_values,
                                    cell_num_nodes,
                                    cell_local_id,
                                    e_group);
      auto grad_phi = GetResidualFFGradPhi(grad_shape_values,
                                           cell_num_nodes,
                                           cell_local_id,
                                           e_group);

      double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi) );

      double rrandom = rng.Rand() * r_abs_cellk_interior_max[cell_local_id];

      //======================================== Determine weight
      if (std::fabs(rrandom) < std::fabs(r))
        new_particle.w = r/std::fabs(r);
      else
      {
        new_particle.w = 0.0;
        new_particle.alive = false;
        particle_rejected = true;
      }

      new_particle.cur_cell_local_id  = cell_local_id;
      new_particle.cur_cell_global_id = cell.global_id;
    }//while particle rejected
  }//interior
    //################################################## SURFACE
  else
  {
    //==================================== Randomly sample face
    int64_t cf = std::lower_bound(surface_cdf.begin(),
                                  surface_cdf.end(),
                                  rng.Rand()) - surface_cdf.begin();

    const auto& rcellface = residual_info_cell_bndry_faces[cf];

    const uint64_t cell_local_id  = rcellface.cell_local_id;
    const auto&    cell           = ref_solver.grid->local_cells[cell_local_id];
    const auto&    cell_pwl_view  = ref_solver.pwl->GetCellMappingFE(cell_local_id);
    const auto&    cell_geom_info = cell_geometry_info->operator[](cell_local_id);
    const size_t   cell_num_nodes = ref_solver.pwl->GetCellNumNodes(cell);

    const int   f    = rcellface.ass_face;
    const auto& face = cell.faces[f];
    const auto& n    = face.normal;

    //======================================== Start rejection sampling
    bool particle_rejected = true;
    while (particle_rejected)
    {
      particle_rejected = false;
      new_particle.alive = true;

      //==================================== Sample position
      new_particle.pos = GetRandomPositionOnCellFace(rng, cell_geom_info, f);

      //==================================== Sample direction

      chi_mesh::Vector3 omega = RandomCosineLawDirection(rng,-1.0*n);
      new_particle.dir = omega;

      //==================================== Populate shape values
      cell_pwl_view->ShapeValues(new_particle.pos, shape_values);

      //==================================== Get Residual
      double phi_P = GetResidualFFPhi(shape_values,
                                      cell_num_nodes,
                                      cell_local_id,
                                      e_group);

      double phi_N = phi_P;
      if (not face.has_neighbor) phi_N = 0.0; //TODO: Specialize for bndries

      double r = (1.0/FOUR_PI)*(phi_N - phi_P);

      double rrandom = rng.Rand()*rcellface.maximum_rstar_absolute;

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