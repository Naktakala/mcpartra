#include "solver_montecarlon.h"

#include "ChiLog/chi_log.h"

extern ChiLog& chi_log;

//###################################################################
/**Initializes raytracers according to cell-size histogram.*/
void mcpartra::Solver::InitRaytracing()
{
  chi_log.Log() << "MCPartra: Initializing Raytracing items.";

  /**Lambda to get cell bounding box.*/
  auto GetCellApproximateSize = [this](const chi_mesh::Cell& cell)
  {
    const auto& v0 = grid->vertices[cell.vertex_ids[0]];
    double xmin = v0.x, xmax = v0.x;
    double ymin = v0.y, ymax = v0.y;
    double zmin = v0.z, zmax = v0.z;

    for (uint64_t vid : cell.vertex_ids)
    {
      const auto& v = grid->vertices[vid];

      xmin = std::min(xmin, v.x); xmax = std::max(xmax, v.x);
      ymin = std::min(ymin, v.y); ymax = std::max(ymax, v.y);
      zmin = std::min(zmin, v.z); zmax = std::max(zmax, v.z);
    }

    return (chi_mesh::Vector3(xmin, ymin, zmin) -
            chi_mesh::Vector3(xmax, ymax, zmax)).Norm();
  };

  cell_sizes.assign(grid->local_cells.size(), 0.0);
  for (const auto& cell : grid->local_cells)
    cell_sizes[cell.local_id] = GetCellApproximateSize(cell);

}