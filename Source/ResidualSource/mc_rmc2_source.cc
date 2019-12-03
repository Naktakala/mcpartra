#include "mc_rmc2_source.h"

#include "../../RandomNumberGenerator/montecarlon_rng.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_polyhedron.h>

#include <PiecewiseLinear/pwl.h>

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
  bool use_uniform_sampling) :
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

  const int ALL_BOUNDRIES = -1;

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

      //Determine if face will be sampled
      bool sample_face = true;


      if (sample_face)
      {
        double face_area = fv_view->face_area[f];
        chi_mesh::Matrix3x3 R;

        chi_mesh::Vector n = cell->faces[f].normal*-1.0;
        chi_mesh::Vector khat(0.0,0.0,1.0);

        if      (n.Dot(khat) >  0.9999999)
          R.SetDiagonalVec(1.0,1.0,1.0);
        else if (n.Dot(khat) < -0.9999999)
          R.SetDiagonalVec(1.0,1.0,-1.0);
        else
        {
          chi_mesh::Vector binorm = khat.Cross(n);
          binorm = binorm/binorm.Norm();

          chi_mesh::Vector tangent = binorm.Cross(n);
          tangent = tangent/tangent.Norm();

          R.SetColJVec(0,tangent);
          R.SetColJVec(1,binorm);
          R.SetColJVec(2,n);
        }

        total_patch_area += face_area;
        source_patches.emplace_back(cell_glob_index,f,R,face_area);
      }
      ++f;
    }//for face
  }//for cell

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
CreateParticle(chi_montecarlon::RandomNumberGenerator* rng)
{
  chi_montecarlon::Particle new_particle;

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
    auto polyh_cell = (chi_mesh::CellPolyhedronV2*)cell;
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
  //============================== Interpolate phi
  auto pwl_view = ref_solver->pwl_discretization->MapFeView(cell_glob_index);
  std::vector<double> shape_values(pwl_view->dofs,0.0);
  pwl_view->ShapeValues(new_particle.pos,shape_values);

  //TODO: Make proper mappings here
  double cell_phi = 0.0;
  for (int dof=0; dof<pwl_view->dofs; dof++)
  {
    int map = (*resid_ff->local_cell_dof_array_address)[cell->cell_local_id];
    int ir = map + dof*resid_ff->num_grps*resid_ff->num_moms +
             resid_ff->num_grps*0 + new_particle.egrp;
    cell_phi += (*resid_ff->field_vector_local)[ir]*shape_values[dof];
  }
  //============================== Get boundary flux
  double bndry_phi = 0.0;
  if (ref_bndry == abs(face.neighbor))
    bndry_phi = 1.0;

  const double FOUR_PI = 4.0*acos(-1.0);

  double omega_dot_n = new_particle.dir.Dot(face.normal*-1.0);

  new_particle.w = -(1.0/FOUR_PI)*(cell_phi - bndry_phi);

  //======================================== Sample energy
  new_particle.egrp = 0;

  new_particle.cur_cell_ind = cell_glob_index;

//  chi_log.Log(LOG_ALL) << new_particle.pos.PrintS();
//  usleep(100000);

  return new_particle;
}
