#include "solver_montecarlon.h"

#include <ChiMesh/Raytrace/raytracing.h>

#include <chi_log.h>
extern ChiLog& chi_log;

//###################################################################
/**The default raytracing algorithm.*/
void chi_montecarlon::Solver::RaytraceRMC(Particle& prtcl)
{
//  //======================================== Get cell
//  chi_mesh::Cell* cell;
//  if (prtcl.cur_cell_local_id >= 0)
//    cell = &(grid->local_cells[prtcl.cur_cell_local_id]);
//  else
//  {
//    cell = grid->cells[prtcl.cur_cell_global_id];
//    prtcl.cur_cell_local_id = cell->local_id;
//  }
//
//  //======================================== Distance to event
//  double d_to_surface = 1.0e15;
//
//  chi_mesh::Vector3 posf = prtcl.pos;
//  chi_mesh::Vector3 dirf = prtcl.dir;
//  int                ef = prtcl.egrp;
//  chi_mesh::RayDestinationInfo ray_dest_info =
//    chi_mesh::RayTrace(grid, cell,
//                       prtcl.pos, prtcl.dir,
//                       d_to_surface, posf);
//
//  //posf set in call to RayTrace
//  ContributeTallyRMC(prtcl,posf,ray_dest_info);
//  if (ray_dest_info.destination_face_neighbor < 0)
//  {
//    bool reflecting = false;
//    if (!reflecting) prtcl.alive = false;
//    else {}
//  }//if bndry
//  else
//  {
//    prtcl.pre_cell_global_id = prtcl.cur_cell_global_id;
//    prtcl.cur_cell_global_id = ray_dest_info.destination_face_neighbor;
//
//    if ((not grid->IsCellLocal(prtcl.cur_cell_global_id)) and
//        (not mesh_is_global))
//    {
//      prtcl.pos = posf;
//      prtcl.dir = dirf;
//      prtcl.egrp = ef;
//      prtcl.banked = true;
//      prtcl.cur_cell_local_id =
//        cell_neighbor_nonlocal_local_id[ray_dest_info.destination_face_neighbor];
//      outbound_particle_bank.push_back(prtcl);
//    }
//    else
//    {
//      int f = ray_dest_info.destination_face_index;
//      int adj_cell_local_id = cell->faces[f].GetNeighborLocalID(grid);
//      prtcl.cur_cell_local_id = adj_cell_local_id;
//    }
//  }
//
//  prtcl.pos = posf;
//  prtcl.dir = dirf;
//  prtcl.egrp = ef;

}
