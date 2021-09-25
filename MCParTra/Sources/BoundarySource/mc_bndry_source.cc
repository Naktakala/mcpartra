#include "mc_bndry_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include <ChiMesh/MeshContinuum/chi_meshcontinuum.h>
#include <ChiMesh/Cell/cell_slab.h>
#include <ChiMesh/Cell/cell_polygon.h>
#include <ChiMesh/Cell/cell_polyhedron.h>

#include <ChiMath/SpatialDiscretization/FiniteVolume/fv.h>
#include <ChiMath/SpatialDiscretization/FiniteVolume/CellViews/fv_slab.h>
#include <ChiMath/SpatialDiscretization/FiniteVolume/CellViews/fv_polygon.h>
#include <ChiMath/SpatialDiscretization/FiniteVolume/CellViews/fv_polyhedron.h>

#include <ChiMath/Statistics/cdfsampler.h>

#include "chi_log.h"

extern ChiLog& chi_log;


//###################################################################
/**Initializes references to necessary mesh elements.
 *
 * Since this is a boundary source the major difficulty here is
 * to get the sampling direction right. Sampled directions need
 * to be rotation to the opposite of the outward pointing boundary
 * normals. For this we assemble a rotation matrix for each face
 * that needs to be sampled.*/
void mcpartra::BoundarySource::
  Initialize(chi_mesh::MeshContinuumPtr&    ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map)
{
  chi_log.Log(LOG_0) << "Initializing Boundary Source";

  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;

  const int ALL_BOUNDRIES = -1;

  //============================================= Build surface src patchess
  // The surface points
  double total_patch_area = 0.0;
  for (auto& cell : grid->local_cells)
  {
    auto fv_view = fv_sdm->MapFeView(cell.local_id);


    int f=0;
    for (auto& face : cell.faces)
    {
      if (face.has_neighbor) {++f; continue;}

      //Determine if face will be sampled
      bool sample_face = false;
      if (ref_bndry == ALL_BOUNDRIES)
        sample_face = true;
      if ((ref_bndry != ALL_BOUNDRIES) and (ref_bndry == face.neighbor_id))
        sample_face = true;

      if (sample_face)
      {
        double face_area = fv_view->face_area[f];
        chi_mesh::Matrix3x3 R;

        chi_mesh::Vector3 n = face.normal*-1.0;
        chi_mesh::Vector3 khat(0.0,0.0,1.0);

        if      (n.Dot(khat) >  0.9999999)
          R.SetDiagonalVec(1.0,1.0,1.0);
        else if (n.Dot(khat) < -0.9999999)
          R.SetDiagonalVec(1.0,1.0,-1.0);
        else
        {
          chi_mesh::Vector3 binorm = khat.Cross(n);
          binorm = binorm/binorm.Norm();

          chi_mesh::Vector3 tangent = binorm.Cross(n);
          tangent = tangent/tangent.Norm();

          R.SetColJVec(0,tangent);
          R.SetColJVec(1,binorm);
          R.SetColJVec(2,n);
        }

        total_patch_area += face_area;
        source_patches.emplace_back(cell.local_id,f,R,face_area);
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
    source_patch_cdf[p] = cumulative_value/total_patch_area;
    ++p;
  }

  int local_num_patches = source_patches.size();
  int global_num_patches = 0;

  MPI_Allreduce(&local_num_patches,
                &global_num_patches,
                1,
                MPI_INT,
                MPI_SUM,
                MPI_COMM_WORLD);

  if (global_num_patches == 0)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::BoundarySource::Initialize"
         " No source patches.";
    exit(EXIT_FAILURE);
  }
}

//###################################################################
/**Source routine for boundary source.*/
mcpartra::Particle mcpartra::BoundarySource::
  CreateParticle(chi_math::RandomNumberGenerator& rng)
{
  mcpartra::Particle new_particle;

  if (source_patch_cdf.empty())
  {
    new_particle.alive = false;
    return new_particle;
  }

  //======================================== Sample source patch
  int source_patch_sample = std::lower_bound(
    source_patch_cdf.begin(),
    source_patch_cdf.end(),
    rng.Rand()) - source_patch_cdf.begin();

  auto& source_patch = source_patches[source_patch_sample];

  //======================================== Get references
  int cell_local_id    = std::get<0>(source_patch);
  int f                = std::get<1>(source_patch);
  auto& RotationMatrix = std::get<2>(source_patch);
  auto& cell           = grid->local_cells[cell_local_id];
  auto& face           = cell.faces[f];
  auto fv_view         = fv_sdm->MapFeView(cell_local_id);

  //======================================== Sample position
  if      (cell.Type() == chi_mesh::CellType::SLAB)
    new_particle.pos = grid->vertices[face.vertex_ids[0]];
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
  {
    chi_mesh::Vertex& v0 = grid->vertices[face.vertex_ids[0]];
    chi_mesh::Vertex& v1 = grid->vertices[face.vertex_ids[1]];
    double w = rng.Rand();
    new_particle.pos = v0*w + v1*(1.0-w);
  }
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
  {
    auto polyh_cell = (chi_mesh::CellPolyhedron*)(&cell);
    auto polyh_fv_view = (PolyhedronFVValues*)fv_view;

    auto edges = polyh_cell->GetFaceEdges(f);

    //===================== Sample side
    double rn = rng.Rand();
    int s = -1;
    double cumulated_area = 0.0;
    for (auto face_side_area : polyh_fv_view->face_side_area[f])
    {
      s++;
      if (rn < ((face_side_area+cumulated_area)/polyh_fv_view->face_area[f]))
        break;
      cumulated_area+=face_side_area;
    }

    double u = rng.Rand();
    double v = rng.Rand();
    while ((u+v)>1.0)
    {u = rng.Rand(); v = rng.Rand();}

    chi_mesh::Vector3& v0 = grid->vertices[edges[s][0]];
    new_particle.pos = v0 + polyh_fv_view->face_side_vectors[f][s][0]*u +
                            polyh_fv_view->face_side_vectors[f][s][1]*v;
  }

  //======================================== Sample direction
  double costheta = rng.Rand();     //Sample half-range only
  double theta    = acos(sqrt(costheta));
  double varphi   = rng.Rand()*2.0*M_PI;

  chi_mesh::Vector3 ref_dir;
  ref_dir.x = sin(theta)*cos(varphi);
  ref_dir.y = sin(theta)*sin(varphi);
  ref_dir.z = cos(theta);

  new_particle.dir = RotationMatrix*ref_dir;

  //======================================== Sample energy
  new_particle.egrp = 0;

  //======================================== Determine weight
  new_particle.w = 1.0;

  new_particle.cur_cell_global_id = cell.global_id;
  new_particle.cur_cell_local_id  = cell.local_id;

  if (ref_solver.options.uncollided_only)
    new_particle.ray_trace_method = mcpartra::RayTraceMethod::UNCOLLIDED;

  return new_particle;
}