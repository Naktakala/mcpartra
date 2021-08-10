#include "raytracing.h"

#include "ChiMesh/Cell/cell.h"
#include "ChiMesh/Cell/cell_polygon.h"
#include "ChiMesh/MeshContinuum/chi_meshcontinuum.h"

#include "chi_log.h"
extern ChiLog& chi_log;


//###################################################################
/**Parent raytracing routine.*/
chi_mesh::RayTracerOutputInformation chi_mesh::RayTracer::
  TraceRay(const Cell &cell,
           Vector3 &pos_i,
           Vector3 &omega_i,
           int function_depth/*=0*/)
{
  RayTracerOutputInformation oi;

  bool intersection_found = false;
  bool backward_tolerance_hit = false;

  if (cell.Type() == chi_mesh::CellType::SLAB)
    TraceSlab(cell, pos_i, omega_i,
              intersection_found    /*byRef*/,
              backward_tolerance_hit/*byRef*/,
              oi                    /*byRef*/);
  else if (cell.Type() == chi_mesh::CellType::POLYGON)
    TracePolygon(cell, pos_i, omega_i,
                 intersection_found    /*byRef*/,
                 backward_tolerance_hit/*byRef*/,
                 oi                    /*byRef*/);
  else if (cell.Type() == chi_mesh::CellType::POLYHEDRON)
    TracePolyhedron(cell, pos_i, omega_i,
                    intersection_found    /*byRef*/,
                    backward_tolerance_hit/*byRef*/,
                    oi                    /*byRef*/);
  else
    throw std::logic_error("Unsupported cell type encountered in call to "
                           "chi_mesh::RayTrace.");

  if (!intersection_found)
  {
    if (function_depth < 5)
    {
      // Nudge particle towards centroid
      chi_mesh::Vector3 v_p_i_cc = (cell.centroid - pos_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      oi = TraceRay(cell,pos_i_nudged,omega_i,function_depth+1);

      return oi;
    }

    if (function_depth < 7)
    {
      // Nudge particle away from line between location and cell center
      chi_mesh::Vector3 v_p_i_cc = (cell.centroid - pos_i).Cross(omega_i);
      chi_mesh::Vector3 pos_i_nudged = pos_i + v_p_i_cc * epsilon_nudge;

      oi = TraceRay(cell,pos_i_nudged,omega_i,function_depth+1);

      return oi;
    }


    std::stringstream outstr;

    outstr
      << "Intersection not found at function level " << function_depth << "."
      << ((backward_tolerance_hit)? " Backward tolerance hit. " : "")
      << "For particle xyz="
      << pos_i.PrintS() << " uvw="
      << omega_i.PrintS() << " " << (pos_i + extension_distance*omega_i).PrintS()
      << " " << extension_distance
      << " in cell " << cell.global_id
      << " with vertices: \n";

    for (auto vi : cell.vertex_ids)
      outstr << grid.vertices[vi].PrintS() << "\n";

    for (auto& face : cell.faces)
    {
      outstr << "Face with centroid: " << face.centroid.PrintS() << " ";
      outstr << "n=" << face.normal.PrintS() << "\n";
      for (auto vi : face.vertex_ids)
        outstr << grid.vertices[vi].PrintS() << "\n";
    }

    outstr << "o Cell\n";
    for (auto& vid : cell.vertex_ids)
    {
      auto& v = grid.vertices[vid];
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (auto& face : cell.faces)
    {
      auto& v = face.centroid;
      outstr << "v " << v.x << " " << v.y << " " << v.z << "\n";
    }

    for (size_t f=0; f < cell.faces.size(); ++f)
    {
      auto& face = cell.faces[f];
      outstr << "f ";
      for (auto vid : face.vertex_ids)
      {
        size_t ref_cell_id=0;
        for (size_t cid=0; cid < cell.vertex_ids.size(); ++cid)
          if (cell.vertex_ids[cid] == vid)
            ref_cell_id = cid + 1;

        outstr << ref_cell_id << "// ";
      }
      outstr << "\n";
    }

    oi.particle_lost = true;
    oi.lost_particle_info = outstr.str();

    std::cout << outstr.str();
    exit(EXIT_FAILURE);
  }

  return oi;
}