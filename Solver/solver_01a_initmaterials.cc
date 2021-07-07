#include"solver_montecarlon.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initialize materials.*/
void mcpartra::Solver::InitMaterials()
{
  chi_log.Log() << "MCParTra: Initializing Materials.";

  typedef chi_physics::TransportCrossSections TrXS;

  //=================================== Check materials exist
  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::Solver::Initialize : No materials found.";
    exit(EXIT_FAILURE);
  }
  size_t num_mat = chi_physics_handler.material_stack.size();

  //=================================== Check that cells all map to a
  //                                    valid material
  size_t invalid_cell_mats = 0;
  for (const auto& cell : grid->local_cells)
    if (cell.material_id < 0 or cell.material_id >= static_cast<int>(num_mat))
      invalid_cell_mats++;

  if (invalid_cell_mats > 0)
    throw std::invalid_argument(
      "MCParTra: " + std::to_string(invalid_cell_mats) + " cells found with "
      "material-IDs pointing to invalid materials.");

  //=================================== Initialize Materials and make property
  //                                    mappings
  matid_xs_map.resize(num_mat,-1);
  matid_q_map.resize(num_mat,-1);
  for (size_t m=0; m<num_mat; m++)
  {
    auto cur_mat = chi_physics_handler.material_stack[m];
    bool material_xs_mapped = false;

    //======================= Only first xs will be used
    size_t num_props = cur_mat->properties.size();
    for (size_t p=0; p<num_props; p++)
    {
      auto& property = cur_mat->properties[p];
      auto property_type = cur_mat->properties[p]->Type();

      if (property_type == chi_physics::PropertyType::TRANSPORT_XSECTIONS and
        (not material_xs_mapped))
      {
        auto transp_xs = std::static_pointer_cast<TrXS>(property);

        transp_xs->ComputeDiscreteScattering(options.scattering_order);

        if (transp_xs->num_groups > num_grps)
          num_grps = transp_xs->num_groups;

        matid_xs_map[m] = static_cast<int>(p);
        material_xs_mapped = true;
      }

      if (property_type == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
        matid_q_map[m] = static_cast<int>(p);
    }//for prop

    using namespace std;

    if (matid_xs_map[m]<0)
      throw invalid_argument(
        "The mapping of material " + to_string(m) + " to a cross section "
        "property failed. This indicates that the given material might not "
        "have a transport cross section property.");
  }//for mat
  MPI_Barrier(MPI_COMM_WORLD);


}