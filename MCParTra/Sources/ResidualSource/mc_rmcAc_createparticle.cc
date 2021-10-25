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
  typedef chi_mesh::Vector3 Vec3;
  constexpr double FOUR_PI = 4.0*M_PI;
  mcpartra::Particle new_particle;

  std::vector<double>            shape_values;
  std::vector<chi_mesh::Vector3> grad_shape_values;

  //======================================== Sample energy group
  size_t g = SampleCDF(group_biased_cdf, rng);
  new_particle.egrp = static_cast<int>(g);

  //======================================== Sample element
  size_t elem = SampleCDF(group_element_biased_cdf[g], rng);

  const auto& src_element = *group_sources[g][elem];
  bool sample_interior = (src_element.type == RessidualInfoType::Interior);

  //======================================== Get reference cell info
  size_t cell_local_id = src_element.cell_local_id;
  const auto&    cell           = ref_solver.grid->local_cells[cell_local_id];
  const auto&    cell_pwl_view  = ref_solver.pwl->GetCellMappingFE(cell_local_id);
  const auto&    cell_geom_info = cell_geometry_info->operator[](cell_local_id);
  const size_t   cell_num_nodes = ref_solver.pwl->GetCellNumNodes(cell);

  const auto imp_info = GetCellImportanceInfo(cell, g);

  const double C_u = src_element.Rstar_absolute;
  const double C_b = src_element.Rstar_psi_star_absolute;

  //################################################## INTERIOR
  if (sample_interior)
  {
    //======================================== Get information
    const auto&  cell_r_info = src_element;

    const VecDbl& nodal_phi = GetResidualFFPhiAtNodes(cell, cell_num_nodes, 0, g);

    //====================================== Get material properties
    MaterialData mat_data;
    PopulateMaterialData(cell.material_id, g, mat_data);

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

      //==================================== Compute Residual
      double phi = GetPhiH(shape_values, nodal_phi, cell_num_nodes);
      Vec3   grad_phi = GetGradPhiH(grad_shape_values, nodal_phi, cell_num_nodes);

      double r = (1.0/FOUR_PI)*( Q - siga*phi - omega.Dot(grad_phi) );
      double pdf_val = r;
      double pdf_random = rng.Rand() * cell_r_info.maximum_rstar_absolute;
      double angular_w_corr = 1.0;

      //==================================== Compute angular bias
//      double psi_star = imp_info.ExpRep(omega);
//      double angular_w_corr = C_b / (C_u * psi_star);
//
//      double pdf_val = r * psi_star;
//      double pdf_random = rng.Rand() * cell_r_info.maximum_rstar_psi_star_absolute;

      //======================================== Determine weight
      if (std::fabs(pdf_random) < std::fabs(pdf_val))
        new_particle.w = (r/std::fabs(r))*
                         group_element_biased_cdf_corr[g][elem]*
                         angular_w_corr;
      else
      {
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
    //====================================== Get information
    const auto& rcellface = static_cast<const RCellFace&>(src_element);

    const VecDbl& nodal_phi = GetResidualFFPhiAtNodes(cell, cell_num_nodes, 0, g);

    const unsigned int f = rcellface.ass_face;
    const auto& face     = cell.faces[f];
    const auto& n        = face.normal;

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
      double phi_P = GetPhiH(shape_values, nodal_phi, cell_num_nodes);

      double phi_N = phi_P;
      if (not face.has_neighbor) phi_N = 0.0; //TODO: Specialize for bndries

      double r = (1.0/FOUR_PI)*(phi_N - phi_P);

      double rrandom = rng.Rand() * rcellface.maximum_rstar_absolute;

      //======================================== Determine weight
      if (std::fabs(rrandom) < std::fabs(r))
        new_particle.w = (r/std::fabs(r))*group_element_biased_cdf_corr[g][elem];
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