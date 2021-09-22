#include "sdsolver.h"

#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include<cmath>

std::pair<double,double> mcpartra::SourceDrivenSolver::
  GetMaterialSigmas(const chi_physics::TransportCrossSections& xs,
                    unsigned int energy_group)
{
  double sigma_t = xs.sigma_t[energy_group];
  double sigma_s = sigma_t - xs.sigma_a[energy_group];

  return {sigma_t, sigma_s};
}

//###################################################################
/**Entry point ray-tracer which selects different ray tracing algorithms.*/
void mcpartra::SourceDrivenSolver::Raytrace(Particle& prtcl)
{
  //======================================== Save prtcl initial info
  const auto pos_i    = prtcl.pos;
  const auto dir_i    = prtcl.dir;
  const auto egrp_i   = prtcl.egrp;
  const auto weight_i = prtcl.w;

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

  //======================================== Raytrace within cell
  default_raytracer->SetTolerancesFromCellSize(cell_sizes[cell->local_id]);
  auto ray_dest_info = default_raytracer->TraceRay(*cell,prtcl.pos,prtcl.dir);

  if (ray_dest_info.particle_lost)
  {
    lost_particles.emplace_back(ray_dest_info.lost_particle_info);
    prtcl.alive = false;
    return;
  }

  double d_to_surface      = ray_dest_info.distance_to_surface;
  const auto& exiting_face = cell->faces[ray_dest_info.destination_face_index];

  //======================================== Get total and scat xs
  auto& xs = *matid_xs_map2[cell->material_id];
  const auto sigma_values = GetMaterialSigmas(xs, egrp_i);

  const double sigma_t = sigma_values.first;
  const double sigma_s = sigma_values.second;

  auto& xs_scattering_cdfs = matid_scattering_cdfs.at(cell->material_id);

  //======================================== Process tally contributions and/or
  //                                         interaction
  bool particle_went_to_surface = false;
  auto pos_f       = pos_i;
  auto dir_f       = dir_i;
  auto egrp_f      = egrp_i;
  auto weight_f    = weight_i;
  bool prtcl_alive = true;

  if (options.uncollided_only)
  {
    particle_went_to_surface = true;
    pos_f = ray_dest_info.pos_f;

    weight_f = ContributeTallyUNC(prtcl, pos_f, sigma_t);
  }

  if (not options.uncollided_only)
  {
    double d_to_intract = -1.0 * log(1.0-rng0.Rand()) / sigma_t;

    bool particle_interacted = (d_to_intract < d_to_surface);
    bool particle_scattered  = (rng0.Rand() < (sigma_s / sigma_t));

    if (particle_interacted and not particle_scattered)
    {
      particle_went_to_surface = false;
      pos_f = pos_i + dir_i * d_to_intract;
      prtcl_alive = false;
    }//if interacts before surface

    if (particle_interacted and particle_scattered)
    {
      particle_went_to_surface = false;
      pos_f = pos_i + dir_i * d_to_intract;

      auto energy_dir = ProcessScattering(prtcl,xs_scattering_cdfs);
      egrp_f = energy_dir.first;
      dir_f  = energy_dir.second;

      if (options.mono_energy && (egrp_f != egrp_i)) prtcl.Kill();
    }//if interacts before surface

    if (not particle_interacted)
    {
      particle_went_to_surface = true;
      pos_f = ray_dest_info.pos_f;
    }

    ContributeTally(prtcl,pos_f);
  }//if not uncollided.only

  if (particle_went_to_surface)
  {
    if (not exiting_face.has_neighbor)
      prtcl_alive = false; //TODO: Add reflecting boundaries

    if (exiting_face.has_neighbor)
    {
      prtcl.pre_cell_global_id = prtcl.cur_cell_global_id;
      prtcl.cur_cell_global_id = ray_dest_info.destination_face_neighbor;

      if (grid->IsCellLocal(prtcl.cur_cell_global_id))
        prtcl.cur_cell_local_id = exiting_face.GetNeighborLocalID(*grid);
    }//if has neighbor
  }//if travelled to surface

  prtcl.pos   = pos_f;
  prtcl.dir   = dir_f;
  prtcl.egrp  = egrp_f;
  prtcl.w     = weight_f;
  prtcl.alive = prtcl_alive;
  if (not local_cell_importance_setting.empty())
    prtcl.pre_cell_importance = prtcl.pre_cell_importance;
}
