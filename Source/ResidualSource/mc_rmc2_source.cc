#include "mc_rmc2_source.h"

#include "../../RandomNumberGenerator/montecarlon_rng.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_polyhedron.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>

#include <ChiMath/Statistics/cdfsampler.h>
#include <ChiMath/Quadratures/quadrature_gausslegendre.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>
#include <MCParTra/Solver/solver_montecarlon.h>

extern ChiLog chi_log;
extern ChiPhysics chi_physics_handler;

#include <tuple>

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualSource2::ResidualSource2(
  chi_physics::FieldFunction *in_resid_ff,
  bool use_uniform_sampling,
  double in_bndry_val) :
  ref_bndry_val(in_bndry_val),
  sample_uniformly(use_uniform_sampling)
{
  type_index = SourceTypes::RESIDUAL;
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
void chi_montecarlon::ResidualSource2::
Initialize(chi_mesh::MeshContinuum *ref_grid,
           SpatialDiscretization_FV *ref_fv_sdm,
           chi_montecarlon::Solver* ref_solver)
{
  chi_log.Log(LOG_0) << "Initializing Residual Bndry Source";
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  this->ref_solver = ref_solver;

  //============================================= Build surface src patchess
  // The surface points
  double total_patch_area = 0.0;
  for (auto& cell_glob_index : grid->local_cell_glob_indices)
  {
    auto cell = grid->cells[cell_glob_index];
    auto fv_view = fv_sdm->MapFeView(cell_glob_index);

    int f=0;
    for (auto& face : cell->faces)
    {
      if (not grid->IsCellBndry(face.neighbor)) {++f; continue;}
      
      double face_area = fv_view->face_area[f];

      chi_mesh::Vector n = cell->faces[f].normal*-1.0;
      chi_mesh::Matrix3x3 R =
        chi_mesh::Matrix3x3::MakeRotationMatrixFromVector(n);

      total_patch_area += face_area;
      source_patches.emplace_back(cell_glob_index,f,R,face_area);

      ++f;
    }//for face
  }//for cell

  BuildCellVolInfo(ref_grid,ref_fv_sdm);

  //============================================= Build source cdf
  source_patch_cdf.resize(source_patches.size(),0.0);
  double cumulative_value = 0.0;
  int p = 0;
  for (auto& source_patch : source_patches)
  {
    cumulative_value += std::get<3>(source_patch);
    source_patch_cdf[p++] = cumulative_value/total_patch_area;
  }
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource2::
  CreateBndryParticle(chi_montecarlon::RandomNumberGenerator* rng)
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
  int cell_glob_index  = std::get<0>(source_patch);
  int f                = std::get<1>(source_patch);
  auto& RotationMatrix = std::get<2>(source_patch);
  auto cell            = grid->cells[cell_glob_index];
  auto face            = cell->faces[f];
  auto fv_view         = fv_sdm->MapFeView(cell_glob_index);

  //======================================== Sample position
  if      (cell->Type() == chi_mesh::CellType::SLAB)
    new_particle.pos = *grid->nodes[face.vertex_ids[0]];
  else if (cell->Type() == chi_mesh::CellType::POLYGON)
  {
    chi_mesh::Vertex& v0 = *grid->nodes[face.vertex_ids[0]];
    chi_mesh::Vertex& v1 = *grid->nodes[face.vertex_ids[1]];
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
    chi_mesh::Vector& v0 = *grid->nodes[edges[s][0]];
    new_particle.pos = v0 + polyh_fv_view->face_side_vectors[f][s][0]*w0 +
                       polyh_fv_view->face_side_vectors[f][s][1]*w1;
  }

  //======================================== Sample direction
  double costheta = rng->Rand();     //Sample half-range only
  double theta    = acos(sqrt(costheta));
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = RotationMatrix*ref_dir;

  //======================================== Determine weight
  // Interpolate phi
  auto pwl_view = ref_solver->pwl_discretization->MapFeView(cell_glob_index);
  std::vector<double> shape_values(pwl_view->dofs,0.0);
  pwl_view->ShapeValues(new_particle.pos,shape_values);

  auto ff_dof_vals = (*resid_ff).GetCellDOFValues(cell->cell_local_id,
                                                  0,new_particle.egrp);

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

  new_particle.cur_cell_ind = cell_glob_index;

  return new_particle;
}

//###################################################################
/**Executes a source sampling for the residual source.*/
chi_montecarlon::Particle chi_montecarlon::ResidualSource2::
CreateParticle(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

  //======================================== Sample cell
//  int lc = std::lower_bound(ref_solver->cell_residual_cdf.begin(),
//                            ref_solver->cell_residual_cdf.end(),
//                            rng->Rand()) -
//                            ref_solver->cell_residual_cdf.begin();

  int num_loc_cells = ref_solver->grid->local_cell_glob_indices.size();
  int lc = std::floor((num_loc_cells)*rng->Rand());

  int cell_glob_index = ref_solver->grid->local_cell_glob_indices[lc];
  auto cell = ref_solver->grid->cells[cell_glob_index];
  auto pwl_view = ref_solver->pwl_discretization->MapFeView(cell_glob_index);
  int map = ref_solver->local_cell_pwl_dof_array_address[lc];

  int mat_id = cell->material_id;
  int xs_id = ref_solver->matid_xs_map[mat_id];

  chi_physics::Material* mat = chi_physics_handler.material_stack[mat_id];
  auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double siga = xs->sigma_ag[0];
  double sigt = xs->sigma_tg[0];
  double sigs = sigt-siga;

  new_particle.cur_cell_ind = cell_glob_index;

  //======================================== Sample position
  new_particle.pos = GetRandomPositionInCell(rng, cell_vol_info[lc]);

  //======================================== Sample direction
  double costheta = 2.0*rng->Rand()-1.0;
  double theta    = acos(costheta);
  double varphi   = rng->Rand()*2.0*M_PI;

  chi_mesh::Vector ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = ref_dir;

  //======================================== Sample uncollided
  std::vector<double> shape_values;
  pwl_view->ShapeValues(new_particle.pos,shape_values);

  double domain_volume = ref_solver->domain_volume;
  double tmf = ref_solver->tally_multipl_factor;

  double weight = 0.0;
  for (int dof=0; dof<pwl_view->dofs; ++dof)
  {
    int ir = map + dof*ref_solver->num_grps*ref_solver->num_moms +
             ref_solver->num_grps*0 + 0;
    weight += shape_values[dof]*ref_solver->phi_pwl_uncollided_rmc[ir];
  }
  weight *= sigs*domain_volume/tmf;

  new_particle.w = weight;
  new_particle.egrp = 0;

  new_particle.alive = true;

  if (std::fabs(weight) < 1.0e-12)
    new_particle.alive = false;

  return new_particle;
}
