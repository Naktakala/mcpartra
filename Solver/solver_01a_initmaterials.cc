#include"solver_montecarlon.h"
#include "ChiPhysics/PhysicsMaterial/transportxsections/material_property_transportxsections.h"

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initialize materials.*/
void chi_montecarlon::Solver::InitMaterials()
{
  chi_log.Log(LOG_0VERBOSE_1) << "Initializing MonteCarlo solver.";

  //=================================== Initialize Materials
  if (chi_physics_handler.material_stack.empty())
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::Solver::Initialize : No materials found.";
    exit(EXIT_FAILURE);
  }
  int num_mat = chi_physics_handler.material_stack.size();
  matid_xs_map.resize(num_mat,-1);
  matid_q_map.resize(num_mat,-1);
  for (int m=0; m<num_mat; m++)
  {
    auto cur_mat = chi_physics_handler.material_stack[m];

    //======================= Only first xs will be used
    size_t num_props = cur_mat->properties.size();
    for (size_t p=0; p<num_props; p++)
    {
      if (cur_mat->properties[p]->Type() ==
          chi_physics::PropertyType::TRANSPORT_XSECTIONS)
      {
        auto transp_xs =
          std::static_pointer_cast<chi_physics::TransportCrossSections>(
            cur_mat->properties[p]);

        transp_xs->ComputeDiffusionParameters();
        transp_xs->ComputeDiscreteScattering(scattering_order);

        if (transp_xs->G > num_grps)
          num_grps = transp_xs->G;

        matid_xs_map[m] = p;
      }

      if (cur_mat->properties[p]->Type() ==
          chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        matid_q_map[m] = p;
      }
    }//for prop
  }//for mat
  MPI_Barrier(MPI_COMM_WORLD);
}