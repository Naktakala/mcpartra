#include "solver_montecarlon.h"

#include <ChiMesh/Raytrace/raytracing.h>

#include <chi_log.h>
extern ChiLog chi_log;

//###################################################################
/**The default raytracing algorithm.*/
void chi_montecarlon::Solver::RaytraceRMC(Particle& prtcl)
{
  //======================================== Get total and scat xs
  auto cell = grid->cells[prtcl.cur_cell_ind];

  //======================================== Distance to event
  double d_to_surface = 1.0e15;

  chi_mesh::Vector posf = prtcl.pos;
  chi_mesh::Vector dirf = prtcl.dir;
  int                ef = prtcl.egrp;
  chi_mesh::RayDestinationInfo ray_dest_info =
    chi_mesh::RayTrace(grid, cell,
                       prtcl.pos, prtcl.dir,
                       d_to_surface, posf);

  //posf set in call to RayTrace
  ContributeTallyRMC(prtcl,posf,ray_dest_info);
  if (ray_dest_info.destination_face_neighbor < 0)
  {
    bool reflecting = false;
    if (!reflecting) prtcl.alive = false;
    else {}
  }//if bndry
  else
  {
    prtcl.pre_cell_ind = prtcl.cur_cell_ind;
    prtcl.cur_cell_ind = ray_dest_info.destination_face_neighbor;
    if ((not mesh_is_global) and (not grid->IsCellLocal(prtcl.cur_cell_ind)))
    {
      prtcl.pos = posf;
      prtcl.dir = dirf;
      prtcl.egrp = ef;
      outbound_particle_bank.push_back(prtcl);
      prtcl.banked = true;
    }
  }

  prtcl.pos = posf;
  prtcl.dir = dirf;
  prtcl.egrp = ef;

}
