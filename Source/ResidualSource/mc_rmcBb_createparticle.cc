#include "mc_rmcB_source.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_polyhedron.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>
#include "../../Solver/solver_montecarlon.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

#include <tuple>


//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSourceB::
CreateBndryParticle(chi_math::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  if (source_patch_cdf.empty())
  {
    new_particle.alive = false;
    return new_particle;
  }

  //======================================== Sample source patch
  int source_patch_sample = std::lower_bound(
    source_patch_cdf.begin(),
    source_patch_cdf.end(),
    rng->Rand()) - source_patch_cdf.begin();

  auto& source_patch = source_patches[source_patch_sample];

  //======================================== Get references
  int cell_local_id    = std::get<0>(source_patch);
  int f                = std::get<1>(source_patch);
  auto& RotationMatrix = std::get<2>(source_patch);
  auto cell            = &grid->local_cells[cell_local_id];
  auto& face           = cell->faces[f];
  auto fv_view         = fv_sdm->MapFeView(cell_local_id);

  //======================================== Sample position
  if      (cell->Type() == chi_mesh::CellType::SLAB)
    new_particle.pos = *grid->vertices[face.vertex_ids[0]];
  else if (cell->Type() == chi_mesh::CellType::POLYGON)
  {
    chi_mesh::Vertex& v0 = *grid->vertices[face.vertex_ids[0]];
    chi_mesh::Vertex& v1 = *grid->vertices[face.vertex_ids[1]];
    double w = rng->Rand();
    new_particle.pos = v0*w + v1*(1.0-w);
  }
  else if (cell->Type() == chi_mesh::CellType::POLYHEDRON)
  {
    auto polyh_cell = (chi_mesh::CellPolyhedron*)cell;
    auto polyh_fv_view = (PolyhedronFVView*)fv_view;

    auto edges = polyh_cell->GetFaceEdges(f);

    //===================== Sample side
    double rn = rng->Rand();
    int s = -1;
    double cumulated_area = 0.0;
    for (auto face_side_area : polyh_fv_view->face_side_area[f])
    {
      s++;
      if (rn < ((face_side_area+cumulated_area)/polyh_fv_view->face_area[f]))
        break;
      cumulated_area+=face_side_area;
    }

    double w0 = rng->Rand();
    double w1 = rng->Rand()*(1.0-w0);
    chi_mesh::Vector3& v0 = *grid->vertices[edges[s][0]];
    new_particle.pos = v0 + polyh_fv_view->face_side_vectors[f][s][0]*w0 +
                       polyh_fv_view->face_side_vectors[f][s][1]*w1;
  }

  //======================================== Sample direction
  double costheta = rng->Rand();     //Sample half-range only
  double theta    = acos(sqrt(costheta));
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = RotationMatrix*ref_dir;

  //======================================== Determine weight
  // Interpolate phi
  auto pwl_view = ref_solver->pwl->MapFeViewL(cell_local_id);
  std::vector<double> shape_values(pwl_view->dofs,0.0);
  pwl_view->ShapeValues(new_particle.pos,shape_values);
  auto ff_dof_vals = (*resid_ff).GetCellDOFValues(cell_local_id,
                                                  new_particle.egrp,
                                                  0);

  double cell_phi = 0.0;
  for (int dof=0; dof<pwl_view->dofs; dof++)
    cell_phi += ff_dof_vals[dof]*shape_values[dof];

  //============================== Get boundary flux
  double bndry_phi = 0.0;
  if (ref_bndry == abs(face.neighbor) and (face.neighbor<0))
    bndry_phi = 2.0*ref_bndry_val;

  new_particle.w = bndry_phi - cell_phi;

  //======================================== Sample energy
  new_particle.egrp = 0;

  new_particle.cur_cell_global_id = cell->global_id;
  new_particle.cur_cell_local_id  = cell_local_id;

  new_particle.ray_trace_method = chi_montecarlon::Solver::RayTraceMethod::UNCOLLIDED;
  new_particle.tally_method = chi_montecarlon::Solver::TallyMethod::RMC_CHAR_RAY;

  return new_particle;
}

//###################################################################
/**Applies the re-execution efforts.*/
bool chi_montecarlon::ResidualSourceB::
  CheckForReExecution()
{
  if (ray_trace_phase)
  {
    ray_trace_phase = false;

    auto& tally_blocks = ref_solver->grid_tally_blocks;
    auto& masks = ref_solver->TallyMaskIndex;
    typedef chi_montecarlon::Solver::TallyMask maskv;

    uncollided_fv_tally = tally_blocks[masks[maskv::DEFAULT_FVTALLY]];
    uncollided_pwl_tally = tally_blocks[masks[maskv::DEFAULT_PWLTALLY]];

    tally_blocks[masks[maskv::DEFAULT_FVTALLY]].ZeroOut();
    tally_blocks[masks[maskv::DEFAULT_PWLTALLY]].ZeroOut();

    ref_solver->DevelopCollidedSource(uncollided_fv_tally);
    return true;
  }
  else
  {
    auto& tally_blocks = ref_solver->grid_tally_blocks;
    auto& masks = ref_solver->TallyMaskIndex;
    typedef chi_montecarlon::Solver::TallyMask maskv;

    tally_blocks[masks[maskv::DEFAULT_FVTALLY]] += uncollided_fv_tally;

    tally_blocks[masks[maskv::DEFAULT_PWLTALLY]] += uncollided_pwl_tally;
    return false;
  }
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSourceB::
CreateCollidedParticle(chi_math::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  //======================================== Sample group
  auto& cdf_group = ref_solver->cdf_phi_unc_group;
  int group_g = std::lower_bound(
    cdf_group.begin(),
    cdf_group.end(),
    rng->Rand()) - cdf_group.begin();

  //======================================== Sample cell
  auto& cdf_cell = ref_solver->cdf_phi_unc_group_cell[group_g];
  int lc = std::lower_bound(
    cdf_cell.begin(),
    cdf_cell.end(),
    rng->Rand()) - cdf_cell.begin();

  auto& tally_blocks = ref_solver->grid_tally_blocks;
  auto& masks = ref_solver->TallyMaskIndex;
  typedef chi_montecarlon::Solver::TallyMask maskv;

  auto& unc_fv_tally  = uncollided_fv_tally.tally_global;
  auto& unc_fem_tally = uncollided_pwl_tally.tally_global;

  //======================================== Randomly Sample Cell
  int cell_glob_index = ref_solver->grid->local_cell_glob_indices[lc];
  auto cell = &ref_solver->grid->local_cells[lc];
  auto pwl_view = ref_solver->pwl->MapFeViewL(lc);
  auto fv_view = ref_solver->fv->MapFeView(lc);

  int mat_id = cell->material_id;
  int xs_id = ref_solver->matid_xs_map[mat_id];

  chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];
  auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double siga = xs->sigma_ag[0];
  double sigt = xs->sigma_tg[0];
  double sigs = sigt-siga;

  //======================================== Sample position
  new_particle.pos = GetRandomPositionInCell(rng, cell_vol_info[lc]);

  //======================================== Sample direction
  double costheta = 2.0*rng->Rand()-1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;

  //======================================== Sample uncollided
  std::vector<double> shape_values;
  pwl_view->ShapeValues(new_particle.pos,shape_values);

  int irfv = ref_solver->fv->MapDOFLocal(cell,&ref_solver->uk_man_fv,/*m*/0,0);
  double weight = 0.0;
  for (int dof=0; dof<pwl_view->dofs; ++dof)
  {
    int irfem = ref_solver->pwl->MapDFEMDOFLocal(cell,dof,&ref_solver->uk_man_fem,/*m*/0,/*g*/0);
    weight += shape_values[dof] * unc_fem_tally[irfem];
  }
  weight *= sigs*ref_solver->IntVSumG_phi_unc/(std::fabs(unc_fv_tally[irfv])/fv_view->volume);

  new_particle.w = weight;
  new_particle.egrp = 0;

  new_particle.alive = true/*true*/;

  if (std::fabs(weight) < 1.0e-12)
    new_particle.alive = false;

  new_particle.cur_cell_global_id = cell_glob_index;
  new_particle.cur_cell_local_id  = cell->local_id;

  new_particle.ray_trace_method = chi_montecarlon::Solver::RayTraceMethod::STANDARD;
  new_particle.tally_method = chi_montecarlon::Solver::TallyMethod::STANDARD;

  return new_particle;
}
