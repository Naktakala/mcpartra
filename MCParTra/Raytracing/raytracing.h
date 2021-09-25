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
  bool         particle_lost = false;
  std::string  lost_particle_info;
};

//###################################################################
/**Raytracer object.*/
class RayTracer
{
private:
  chi_mesh::MeshContinuumPtr reference_grid;
  std::vector<double> cell_sizes;
public:
  double epsilon_nudge      = 1.0e-8;
  double backward_tolerance = 1.0e-10;
  double extension_distance = 1.0e5;
  bool   perform_concavity_checks = true;

  explicit
  RayTracer(double in_epsilon_nudge      = 1.0e-8,
            double in_backward_tolerance = 1.0e-10,
            double in_extension_distance = 1.0e5,
            bool   in_perform_concavity_checks = true) :
    epsilon_nudge     (in_epsilon_nudge     ),
    backward_tolerance(in_backward_tolerance),
    extension_distance(in_extension_distance),
    perform_concavity_checks(in_perform_concavity_checks)
  {}

private:
  const chi_mesh::MeshContinuum& Grid() const;

  void SetTolerancesFromCellSize(double cell_size)
  {
    epsilon_nudge = cell_size * 1.0e-2;
    backward_tolerance = cell_size * 1.0e-10;
    extension_distance = 3.0 * cell_size;
  }

public:
  void SetGrid(chi_mesh::MeshContinuumPtr& input_grid)
  {reference_grid = input_grid;}

  void SetCellSizes(const std::vector<double>& input)
  {cell_sizes = input;}

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