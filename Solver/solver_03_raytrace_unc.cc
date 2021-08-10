#include "solver_montecarlon.h"

#include "ChiMesh/Raytrace/raytracing.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include<cmath>

//###################################################################
/**The default raytracing algorithm.*/
void mcpartra::Solver::RaytraceUNC(Particle& prtcl)
{
  //======================================== Get cell
  chi_mesh::Cell* cell;
  if (prtcl.cur_cell_local_id >= 0)
  {
    cell = &(grid->local_cells[prtcl.cur_cell_local_id]);
    prtcl.cur_cell_global_id = cell->global_id;
  }
  else
  {
    cell = &grid->cells[prtcl.cur_cell_global_id];
    prtcl.cur_cell_local_id = cell->local_id;
    prtcl.cur_cell_global_id = cell->global_id;
  }

  //======================================== Process cell importance
  if (not local_cell_importance_setting.empty())
  {
    if (prtcl.cur_cell_global_id != prtcl.pre_cell_global_id)
    {
      prtcl.cur_cell_importance = local_cell_importance[cell->local_id];
      ProcessImportanceChange(prtcl);
      if ((prtcl.banked) or (not prtcl.alive)) return;
    }
  }

  //======================================== Get total xs
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  auto mat = chi_physics_handler.material_stack[mat_id];

  auto xs = std::static_pointer_cast<chi_physics::TransportCrossSections>(
    mat->properties[xs_id]);

  double sigt = xs->sigma_t[prtcl.egrp];

  //======================================== Compute distance to event
  double d_to_surface = 1.0e15;

  chi_mesh::Vector3 posf = prtcl.pos;
  chi_mesh::Vector3 dirf = prtcl.dir;
  int                ef = prtcl.egrp;
//  chi_mesh::RayDestinationInfo ray_dest_info =
//    chi_mesh::RayTrace(*grid,          //[Input] Grid
//                       *cell,          //[Input] Current cell
//                       prtcl.pos,     //[Input] Current position
//                       prtcl.dir,     //[Input] Current direction
//                       d_to_surface,  //[Otput] Distance to next surface
//                       posf);         //[Otput] Intersection point at next surf

  auto ray_dest_info = default_raytracer->TraceRay(*cell,prtcl.pos,prtcl.dir);

  if (ray_dest_info.particle_lost)
  {
    lost_particles.emplace_back(ray_dest_info.lost_particle_info);
    prtcl.alive = false;
    return;
  }

  d_to_surface = ray_dest_info.distance_to_surface;
  posf = ray_dest_info.pos_f;

  //======================================== Process surface
  if (d_to_surface <0.0)
  {
    chi_log.Log(LOG_ALLERROR) << "Negative distance to surface.";
    exit(EXIT_FAILURE);
  }

  //posf set in call to RayTrace
  if (prtcl.tally_method == TallyMethod::STANDARD)
    ContributeTallyUNC(prtcl,posf,sigt);
//  if (prtcl.tally_method == TallyMethod::RMC_CHAR_RAY)
//    ContributeTallyRMC(prtcl,posf,ray_dest_info);

  //======================= If surface is boundary
  if (not cell->faces[ray_dest_info.destination_face_index].has_neighbor)
  {
    //TODO: Begin - Add reflecting boundaries
    bool reflecting = false;
    if (!reflecting) prtcl.alive = false;
    else {}
    //TODO: End - Add reflecting boundaries
  }//if bndry
    //======================= If surface is cell face
  else
  {
    prtcl.pre_cell_global_id = prtcl.cur_cell_global_id;
    prtcl.cur_cell_global_id = ray_dest_info.destination_face_neighbor;

    if (not grid->IsCellLocal(prtcl.cur_cell_global_id))
    {
      prtcl.pos = posf;
      prtcl.dir = dirf;
      prtcl.egrp = ef;
      prtcl.banked = true;
      prtcl.cur_cell_local_id =
        cell_neighbor_nonlocal_local_id[ray_dest_info.destination_face_neighbor];
      outbound_particle_bank.push_back(prtcl);
    }
    else
    {
      int f = ray_dest_info.destination_face_index;
      int adj_cell_local_id = cell->faces[f].GetNeighborLocalID(*grid);
      prtcl.cur_cell_local_id = adj_cell_local_id;
    }
  }//not to bndry

  prtcl.pos = posf;
  prtcl.dir = dirf;
  prtcl.egrp = ef;
  if (not local_cell_importance_setting.empty())
    prtcl.pre_cell_importance = prtcl.pre_cell_importance;
}
