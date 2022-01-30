#include"sdsolver.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include <assert.h>

//###################################################################
/**Initialize materials.*/
void mcpartra::SourceDrivenSolver::InitMaterials()
{
  chi_log.Log() << "MCParTra: Initializing Materials.";

  typedef chi_physics::TransportCrossSections TrXS;
  typedef chi_physics::IsotropicMultiGrpSource IsoMGSrc;

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
  matid_has_q_flags.assign(num_mat, false);
  for (size_t m=0; m<num_mat; m++)
  {
    auto cur_mat = chi_physics_handler.material_stack[m];
    bool material_xs_mapped = false;
    int mat_id = static_cast<int>(m);

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

        if (transp_xs->num_groups > num_groups)
          num_groups = transp_xs->num_groups;

        matid_xs_map2[mat_id] = transp_xs;

        auto L = options.scattering_order;
        matid_scattering_cdfs.insert(
          std::make_pair(m,MultigroupScatteringCDFs(*transp_xs, L+1)));
        material_xs_mapped = true;
      }

      if (property_type == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        auto iso_mg_src = std::static_pointer_cast<IsoMGSrc>(property);
        matid_has_q_flags[mat_id] = true;
        matid_q_map2[mat_id]      = iso_mg_src;
      }
    }//for prop

    using namespace std;

    if (not material_xs_mapped)
      throw invalid_argument(
        "The mapping of material " + to_string(m) + " to a cross section "
        "property failed. This indicates that the given material might not "
        "have a transport cross section property.");
  }//for mat

  chi_log.Log() << "MCParTra: Number of groups = " << num_groups;
  MPI_Barrier(MPI_COMM_WORLD);

//  exit(EXIT_SUCCESS);
}


