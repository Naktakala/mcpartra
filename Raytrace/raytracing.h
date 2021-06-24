#ifndef CHI_MESH_RAYTRACING2_H
#define CHI_MESH_RAYTRACING2_H

#include "ChiMesh/Raytrace/raytracing.h"

namespace chi_mesh
{
struct RayTracerOutputInformation
{
  double       distance_to_surface = 0.0;
  Vector3      pos_f;
  unsigned int destination_face_index = 0;
  uint64_t     destination_face_neighbor = 0;
};

//###################################################################
/**Raytracer object.*/
class RayTracer
{
private:
  const chi_mesh::MeshContinuum& grid;
public:
  double epsilon_nudge      = 1.0e-8;
  double backward_tolerance = 1.0e-10;
  double extension_distance = 1.0e5;
  bool   perform_concavity_checks = true;

  explicit
  RayTracer(const chi_mesh::MeshContinuum& in_grid,
            double in_epsilon_nudge      = 1.0e-8,
            double in_backward_tolerance = 1.0e-10,
            double in_extension_distance = 1.0e5,
            bool   in_perform_concavity_checks = false) :
    grid(in_grid),
    epsilon_nudge     (in_epsilon_nudge     ),
    backward_tolerance(in_backward_tolerance),
    extension_distance(in_extension_distance),
    perform_concavity_checks(in_perform_concavity_checks)
  {}

  RayTracerOutputInformation
    TraceRay(const Cell& cell,
             Vector3& pos_i,
             Vector3& omega_i,
             int function_depth=0);

private:
  void TraceSlab(const Cell& cell,
                 Vector3& pos_i,
                 Vector3& omega_i,
                 bool& intersection_found,
                 bool& backward_tolerance_hit,
                 RayTracerOutputInformation& oi);
  void TracePolygon(const Cell& cell,
                    Vector3& pos_i,
                    Vector3& omega_i,
                    bool& intersection_found,
                    bool& backward_tolerance_hit,
                    RayTracerOutputInformation& oi);
  void TracePolyhedron(const Cell& cell,
                       Vector3& pos_i,
                       Vector3& omega_i,
                       bool& intersection_found,
                       bool& backward_tolerance_hit,
                       RayTracerOutputInformation& oi);
};
}

#endif //CHI_MESH_RAYTRACING2_H