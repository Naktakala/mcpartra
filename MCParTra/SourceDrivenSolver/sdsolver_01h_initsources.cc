#include "sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Cell geometry data comprises lines, triangles and tetrahedra that
 * are essential to break cells into their most basic constituents.
 * This routine builds and stores this data.*/
void mcpartra::SourceDrivenSolver::InitCellGeometryData()
{
  chi_log.Log() << "MCParTra: Initializing cell geometry data";
  size_t num_local_cells = grid->local_cells.size();

  cell_geometry_info.reserve(num_local_cells);
  for (const auto& cell : grid->local_cells)
  {
    CellGeometryData cell_info;

    //====================================== Make volume elements
    cell_info.volume_elements = GetCellVolumeElements(cell, grid);
    for (const auto& element : cell_info.volume_elements)
      cell_info.total_volume += element.Volume();

    //====================================== Build volume_element cdf
    {
      const double V = cell_info.total_volume;
      auto& cdf = cell_info.volume_elements_cdf;
      const size_t num_elems = cell_info.volume_elements.size();

      cdf.reserve(num_elems);

      double incr_volume = 0.0;
      for (const auto& element : cell_info.volume_elements)
      {incr_volume += element.Volume(); cdf.push_back(incr_volume / V);}
    }

    //====================================== Make surface elements
    size_t num_faces = cell.faces.size();
    cell_info.faces_surface_elements.reserve(num_faces);
    cell_info.faces_total_area.reserve(num_faces);
    double total_surface_area = 0.0;
    for (const auto& face : cell.faces)
    {
      auto face_surf_elements = GetCellSurfaceSourceElements(cell, face, grid);

      double face_area = 0.0;
      for (const auto& element : face_surf_elements)
        face_area += element.Area();

      total_surface_area += face_area;

      cell_info.faces_total_area.push_back(face_area);
      cell_info.faces_surface_elements.push_back(std::move(face_surf_elements));
    }//for face
    cell_info.total_area = total_surface_area;

    //====================================== Build faces cdf
    {
      const auto& areas = cell_info.faces_total_area;
      auto& cdf = cell_info.faces_cdf;
      const double total_area = cell_info.total_area;

      cdf.reserve(num_faces);

      double incr_area = 0.0;
      for (const auto& area : areas)
      {incr_area += area; cdf.push_back(incr_area / total_area);}
    }

    //====================================== Build surface element cdfs
    cell_info.faces_surface_elements_cdf.resize(num_faces);
    for (size_t f=0; f<num_faces; ++f)
    {
      const auto& areas = cell_info.faces_total_area;
      const auto& elements = cell_info.faces_surface_elements[f];
      auto& cdf = cell_info.faces_surface_elements_cdf[f];
      size_t num_elems = elements.size();

      cdf.reserve(num_elems);

      double incr_area = 0.0;
      for (const auto& element : elements)
      {incr_area += element.Area(); cdf.push_back(incr_area / areas[f]);}
    }//for face f

    cell_geometry_info.push_back(cell_info);
  }//for local cell

}

//###################################################################
/**Initialize Sources.*/
void mcpartra::SourceDrivenSolver::InitSources()
{
  chi_log.Log() << "MCParTra: Initializing sources";

  //============================================= Initialize importance
  size_t num_local_cells = grid->local_cells.size();
  size_t num_imp_elements = num_local_cells *  num_groups;
  local_cell_importance_info.resize(num_imp_elements);
  for (auto& info : local_cell_importance_info) info.importance = 1.0;

  //============================================= Map unverified importance info
  const auto& importance_map = unverified_local_cell_importance_info;
  for (auto& cg_index_info_pair : importance_map)
  {
    auto& cg_index_pair = cg_index_info_pair.first;
    auto& info          = cg_index_info_pair.second;

    uint64_t cell_global_id = cg_index_pair.first;
    uint     group_index    = cg_index_pair.second;

    if (grid->IsCellLocal(cell_global_id))
    {
      const auto& cell = grid->cells[cell_global_id];

      size_t mg_map = cell.local_id * num_groups + group_index;

      local_cell_importance_info[mg_map] = info;
    }//if local
  }//for record

  //============================================= Initialize sources
  for (auto& source : sources)
  {
    source->Initialize(grid,
                       fv,
                       num_groups,
                       m_to_ell_em_map,
                       cell_geometry_info);

    total_globl_source_rate += source->GlobalSourceRate();
    total_local_source_rate += source->LocalSourceRate();
  }
  source_normalization = total_globl_source_rate;

  total_globl_source_rate = 0.0;
  total_local_source_rate = 0.0;
  for (auto& source : sources)
  {
    source->BiasCDFs(options.apply_source_importance_sampling);
    total_globl_source_rate += source->GlobalSourceRate();
    total_local_source_rate += source->LocalSourceRate();
  }
  source_normalization = total_globl_source_rate;


  char buffer[100];
  snprintf(buffer,100,"%g",total_globl_source_rate);
  std::string stotal_globl_src_rate(buffer);

  chi_log.Log() << "  Total amount of sources     = " << sources.size();
  chi_log.Log() << "  Total source rate           = " << stotal_globl_src_rate;

  //============================================= Build local cdf
  size_t num_sources = sources.size();
  local_source_cdf.assign(num_sources, 0.0);
  double running_total = 0.0;
  for (size_t s=0; s<num_sources; ++s)
  {
    running_total += sources[s]->LocalSourceRate();
    local_source_cdf[s] = running_total/total_local_source_rate;
  }

}