#include"sdsolver.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;


//###################################################################
/**Initialize field functions.*/
void mcpartra::SourceDrivenSolver::InitFieldFunctions()
{
  chi_log.Log() << "MCParTra: Initializing Field Functions";
  typedef SpatialDiscretization SD;

  auto fv_sd = std::dynamic_pointer_cast<SD>(fv);
  for (size_t g=0; g < num_groups; g++)
  {
    for (size_t m=0; m < num_moments; m++)
    {
      std::string text_name = TextName() +
                              std::string("-FVFlux_g") +
                              std::to_string(g) +
                              std::string("_m") + std::to_string(m);

      auto group_ff = std::make_shared<chi_physics::FieldFunction>(
        text_name,                                     //Text name
        fv_sd,                                         //Spatial Discretization
        &grid_tally_blocks[TallyMaskIndex[DEFAULT_FVTALLY]].tally_global,
        uk_man_fv,                                     //Unknown manager
        m, g);                                         //Reference unknown and component

      chi_physics_handler.fieldfunc_stack.push_back(group_ff);
      field_functions.push_back(group_ff);
    }//for m
  }//for g

//  if (make_pwld)
  {
    auto pwl_sd = std::dynamic_pointer_cast<SD>(pwl);
    for (size_t g=0; g < num_groups; g++)
    {
      for (size_t m=0; m < num_moments; m++)
      {
        std::string text_name = TextName() +
                                std::string("-PWLFlux_g") +
                                std::to_string(g) +
                                std::string("_m") + std::to_string(m);

        auto group_ff = std::make_shared<chi_physics::FieldFunction>(
          text_name,                                    //Text name
          pwl_sd,                                       //Spatial Discretization
          &grid_tally_blocks[TallyMaskIndex[DEFAULT_PWLTALLY]].tally_global,
          uk_man_pwld,                            //Nodal variable structure
          m, g);                                         //Reference variable and component

        chi_physics_handler.fieldfunc_stack.push_back(group_ff);
        field_functions.push_back(group_ff);
      }//for m
    }//for g
  }//if make_pwld
}