#include"solver_montecarlon.h"
#include <ChiPhysics/chi_physics.h>

extern ChiPhysics&  chi_physics_handler;

//###################################################################
/**Initialize field functions.*/
void chi_montecarlon::Solver::InitFieldFunctions()
{
  for (int g=0; g<num_grps; g++)
  {
    std::string text_name = std::string("Flux_g") + std::to_string(g);

    auto group_ff = new chi_physics::FieldFunction(
      text_name,                                    //Text name
      fv,                                           //Spatial Discretization
      &grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].tally_sigma,
      uk_man_fv,                                    //Unknown manager
      0,g);                                         //Reference unknown and component


//    auto group_ff = new chi_physics::FieldFunction(
//      text_name,
//      chi_physics::FieldFunctionType::FV,
//      fv,
//      &uk_man_fv,
//      0,
//      &grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].tally_global
//      );

    chi_physics_handler.fieldfunc_stack.push_back(group_ff);
    field_functions.push_back(group_ff);
  }

//  if (make_pwld)
  {
    for (int g=0; g<num_grps; g++)
    {
      for (int m=0; m<num_moms; m++)
      {
        std::string text_name = std::string("Flux_g") +
                                std::to_string(g) +
                                std::string("_m") + std::to_string(m);

        auto group_ff = new chi_physics::FieldFunction(
          text_name,                                    //Text name
          chi_physics_handler.fieldfunc_stack.size(),   //FF-id
          chi_physics::FieldFunctionType::DFEM_PWL,     //Type
          grid,                                         //Grid
          pwl,                           //Spatial Discretization
          num_grps,                                     //Number of components
          num_moms,                                     //Number of sets
          g, m,                                          //Ref component, ref set
          &pwl->cell_dfem_block_address,            //Dof block address
//          &phi_pwl_global                           //Data vector
          &grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].tally_global
          );

//        auto group_ff = new chi_physics::FieldFunction(
//          text_name,
//          chi_physics::FieldFunctionType::DFEM_PWL,
//          pwl,
//          &uk_man_fem,
//          0,
//          &grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].tally_global
//          );

        chi_physics_handler.fieldfunc_stack.push_back(group_ff);
        field_functions.push_back(group_ff);
      }//for m
    }//for g
  }//if make_pwld
}