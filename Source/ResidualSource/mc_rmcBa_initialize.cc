#include "mc_rmcB_source.h"

#include <ChiMesh/Cell/cell_polyhedron.h>

#include <FiniteVolume/fv.h>
#include <FiniteVolume/CellViews/fv_polyhedron.h>

#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>
#include <ChiPhysics/PhysicsMaterial/property11_isotropic_mg_src.h>

#include <ChiMesh/FieldFunctionInterpolation/chi_ffinterpolation.h>

#include <ChiPhysics/chi_physics.h>
#include <chi_log.h>
#include "../../Solver/solver_montecarlon.h"

extern ChiLog& chi_log;
extern ChiPhysics&  chi_physics_handler;

#include <tuple>

//###################################################################
/**Constructor for residual source.*/
chi_montecarlon::ResidualSourceB::
ResidualSourceB(chi_physics::FieldFunction *in_resid_ff,
                bool use_uniform_sampling,
                double in_bndry_val) :
  ref_bndry_val(in_bndry_val),
  sample_uniformly(use_uniform_sampling)
{
  type_index = SourceTypes::RESIDUAL_TYPE_B;
  resid_ff = in_resid_ff;
}

//###################################################################
/**Initializes an rmc source.
 *
 * This process involves numerous steps. One of the first steps is
 * to */
void chi_montecarlon::ResidualSourceB::
Initialize(chi_mesh::MeshContinuum *ref_grid,
           SpatialDiscretization_FV *ref_fv_sdm,
           chi_montecarlon::Solver* ref_solver)
{
  chi_log.Log(LOG_0) << "Initializing Residual Bndry Source";
  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  this->ref_solver = ref_solver;

  //============================================= Build surface src patches
  // The surface points
  double total_patch_area = 0.0;
  for (auto& cell : grid->local_cells)
  {
    auto fv_view = fv_sdm->MapFeView(cell.local_id);

    int f=0;
    for (auto& face : cell.faces)
    {
      if (not grid->IsCellBndry(face.neighbor)) {++f; continue;}

      double face_area = fv_view->face_area[f];

      chi_mesh::Vector3 n = face.normal*-1.0;
      chi_mesh::Matrix3x3 R =
        chi_mesh::Matrix3x3::MakeRotationMatrixFromVector(n);

      total_patch_area += face_area;
      source_patches.emplace_back(cell.local_id,f,R,face_area);

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
/**Gets the relative source strength accross all processors.*/
double chi_montecarlon::ResidualSourceB::GetRMCParallelRelativeSourceWeight()
{
  if (ray_trace_phase)
  {
    double local_surface_area = 0.0;
    for (auto& source_patch : source_patches)
      local_surface_area += std::get<3>(source_patch);

    double global_surface_area = 0.0;
    MPI_Allreduce(&local_surface_area,  //sendbuf
                  &global_surface_area, //recvbuf
                  1,                    //recvcount
                  MPI_DOUBLE,           //datatype
                  MPI_SUM,              //operation
                  MPI_COMM_WORLD);      //communicator

    return local_surface_area/global_surface_area;
  }
  else
    return 1.0;
}
