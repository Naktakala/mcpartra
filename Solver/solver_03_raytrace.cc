#include "solver_montecarlon.h"

#include "ChiMesh/Raytrace/raytracing.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include<cmath>

//###################################################################
/**Entry point ray-tracer which selects different ray tracing algorithms.*/
void chi_montecarlon::Solver::Raytrace(Particle& prtcl)
{
  switch (prtcl.ray_trace_method)
  {
    case RayTraceMethod::STANDARD:
      RaytraceSTD(prtcl);
      break;

    case RayTraceMethod::UNCOLLIDED:
      RaytraceUNC(prtcl);
      break;

    default:
      RaytraceSTD(prtcl);
      break;
  }
}

//###################################################################
/**Standard neutral particle raytracer, tracing to next event
 * or to next surface.*/
void chi_montecarlon::Solver::RaytraceSTD(Particle& prtcl)
{
//  chi_log.Log() << "here " << prtcl.pos.PrintS() << " " << prtcl.dir.PrintS() << " " << prtcl.cur_cell_global_id;
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

  //======================================== Get total and scat xs
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  auto mat = chi_physics_handler.material_stack[mat_id];

  auto xs =
    std::static_pointer_cast<chi_physics::TransportCrossSections>(
      mat->properties[xs_id]);

  double sigt = xs->sigma_tg[prtcl.egrp];
  double sigs = sigt - xs->sigma_ag[prtcl.egrp];

  //======================================== Compute distance to event
  double d_to_intract = -1.0*log(1.0-rng0.Rand())/sigt;
  double d_to_surface = 1.0e15;

  chi_mesh::Vector3 posf = prtcl.pos;
  chi_mesh::Vector3 dirf = prtcl.dir;
  int                ef = prtcl.egrp;
  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(*grid,          //[Input] Grid
                       *cell,          //[Input] Current cell
                       prtcl.pos,     //[Input] Current position
                       prtcl.dir,     //[Input] Current direction
                       d_to_surface,  //[Otput] Distance to next surface
                       posf);         //[Otput] Intersection point at next surf

  //======================================== Process interaction
  if (d_to_intract < d_to_surface)
  {
    posf = prtcl.pos + prtcl.dir*d_to_intract;

    if (uncollided_only)
    {
      ef = prtcl.egrp;
      dirf = prtcl.dir;
      prtcl.w *= ((sigs/sigt));
      prtcl.alive = true;
    }
    else
    {
      if (rng0.Rand() < (sigs/sigt))
      {
        auto energy_dir = ProcessScattering(prtcl,xs);
        ef   = energy_dir.first;
        dirf = energy_dir.second;

        if (mono_energy && (ef != prtcl.egrp))
          prtcl.alive = false;
      }
      else
        prtcl.alive = false;

    }

    if (prtcl.tally_method == TallyMethod::STANDARD)
      ContributeTally(prtcl,posf);
  }
    //======================================== Process surface
  else
  {
    if (d_to_surface <0.0)
    {
      chi_log.Log(LOG_ALLERROR) << "Negative distance to surface.";
      exit(EXIT_FAILURE);
    }

    //posf set in call to RayTrace
    if (prtcl.tally_method == TallyMethod::STANDARD)
      ContributeTally(prtcl,posf);

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

      if ((not mesh_is_global) and (not grid->IsCellLocal(prtcl.cur_cell_global_id)))
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
  }//trace to surface
  prtcl.pos = posf;
  prtcl.dir = dirf;
  prtcl.egrp = ef;
  if (not local_cell_importance_setting.empty())
    prtcl.pre_cell_importance = prtcl.pre_cell_importance;
}


