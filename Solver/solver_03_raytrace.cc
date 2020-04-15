#include "solver_montecarlon.h"

#include <ChiMesh/Raytrace/raytracing.h>
#include <ChiPhysics/PhysicsMaterial/property10_transportxsections.h>

#include <chi_log.h>
extern ChiLog chi_log;

#include <ChiPhysics/chi_physics.h>
extern ChiPhysics chi_physics_handler;

#include<cmath>

//###################################################################
/**The default raytracing algorithm.*/
void chi_montecarlon::Solver::Raytrace(Particle& prtcl)
{
  //======================================== Get cell
  chi_mesh::Cell* cell;
  if (prtcl.cur_cell_local_id >= 0)
    cell = &(grid->local_cells[prtcl.cur_cell_local_id]);
  else
  {
    cell = grid->cells[prtcl.cur_cell_global_id];
    prtcl.cur_cell_local_id = cell->local_id;
  }

  //======================================== Get total and scat xs
  int mat_id = cell->material_id;
  int xs_id = matid_xs_map[mat_id];

  auto mat = chi_physics_handler.material_stack[mat_id];

  auto xs = (chi_physics::TransportCrossSections*)mat->properties[xs_id];

  double sigt = xs->sigma_tg[prtcl.egrp];
  double sigs = sigt - xs->sigma_ag[prtcl.egrp];

  //======================================== Compute distance to event
  double d_to_intract = -1.0*log(1.0-rng0.Rand())/sigt;
  double d_to_surface = 1.0e15;

  chi_mesh::Vector3 posf = prtcl.pos;
  chi_mesh::Vector3 dirf = prtcl.dir;
  int                ef = prtcl.egrp;
  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(grid, cell,
                       prtcl.pos, prtcl.dir,
                       d_to_surface, posf);

  //======================================== Process interaction
  if (d_to_intract < d_to_surface)
  {
    posf = prtcl.pos + prtcl.dir*d_to_intract;

    if (!uncollided_only)
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
    else
    {
      ef = prtcl.egrp;
      dirf = prtcl.dir;
      prtcl.w *= ((sigs/sigt));
      prtcl.alive = true;
    }

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
    ContributeTally(prtcl,posf);
    if (ray_dest_info.destination_face_neighbor < 0)
    {
      bool reflecting = false;
      if (!reflecting) prtcl.alive = false;
      else {}
    }//if bndry
    else
    {
      prtcl.cur_cell_global_id = ray_dest_info.destination_face_neighbor;

      if ((not mesh_is_global) and (not grid->IsCellLocal(prtcl.cur_cell_global_id)))
      {
        prtcl.pos = posf;
        prtcl.dir = dirf;
        prtcl.egrp = ef;
        prtcl.banked = true;
//        prtcl.cur_cell_local_id = -1;
        prtcl.cur_cell_local_id =
          cell_neighbor_nonlocal_local_id[ray_dest_info.destination_face_neighbor];
        outbound_particle_bank.push_back(prtcl);
      }
      else
      {
        int f = ray_dest_info.destination_face_index;
        int adj_cell_local_id = cell->faces[f].GetNeighborLocalID(grid);
        prtcl.cur_cell_local_id = adj_cell_local_id;
      }
    }
  }
  prtcl.pos = posf;
  prtcl.dir = dirf;
  prtcl.egrp = ef;

}
