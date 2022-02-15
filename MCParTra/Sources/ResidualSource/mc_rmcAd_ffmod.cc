#include "mc_rmcA_source.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

void mcpartra::ResidualSourceA::RemoveFFDiscontinuities()
{
  typedef SpatialDiscretization_PWLD PWLD;
  //======================================== Check correct ff type
  if (resid_ff->spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::ResidualSource3: "
      << "Only PWLD spatial discretizations are supported for now.";
    exit(EXIT_FAILURE);
  }

  const auto& pwl =
    std::dynamic_pointer_cast<PWLD>(resid_ff->spatial_discretization);
  auto& uk_man = resid_ff->unknown_manager;

  auto pwl_cfem = SpatialDiscretization_PWLC::New(grid,
                    chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  auto cfem_num_local_dofs = pwl_cfem->GetNumLocalDOFs(uk_man);
  auto cfem_num_globl_dofs = pwl_cfem->GetNumGlobalDOFs(uk_man);

  Vec x = chi_math::PETScUtils::CreateVector(
    static_cast<int64_t>(cfem_num_local_dofs),
    static_cast<int64_t>(cfem_num_globl_dofs));

  VecSet(x,0.0);
  Vec x_count;
  VecDuplicate(x,&x_count);

  std::vector<int64_t> global_cfem_dofs;
  for (auto& cell : grid->local_cells)
  {
    for (size_t i=0; i<cell.vertex_ids.size(); ++i)
      for (size_t m=0; m<uk_man.unknowns.size(); ++m)
        for (size_t g=0; g<uk_man.unknowns[m].num_components; ++g)
        {
          int64_t ir_dfem = pwl->MapDOFLocal(cell,i,uk_man,m,g);
          int64_t ir_cfem = pwl_cfem->MapDOF(cell,i,uk_man,m,g);

          global_cfem_dofs.push_back(ir_cfem);

          double val = (*resid_ff->field_vector_local)[ir_dfem];

          VecSetValue(x      ,ir_cfem,val,ADD_VALUES);
          VecSetValue(x_count,ir_cfem,1.0,ADD_VALUES);
        }
  }

  VecAssemblyBegin(x);
  VecAssemblyBegin(x_count);
  VecAssemblyEnd(x);
  VecAssemblyEnd(x_count);

  VecPointwiseDivide(x,x,x_count); // x[i] = x[i]/x_count[i]

  std::vector<double> local_data;
  chi_math::PETScUtils::CopyGlobalVecToSTLvector(x,global_cfem_dofs,local_data);

  uint64_t ic=0;
  for (auto& cell : grid->local_cells)
  {
    for (size_t i=0; i<cell.vertex_ids.size(); ++i)
      for (size_t m=0; m<uk_man.unknowns.size(); ++m)
        for (size_t g=0; g<uk_man.unknowns[m].num_components; ++g)
        {
          int64_t ir_dfem = pwl->MapDOFLocal(cell,i,uk_man,m,g);

          (*resid_ff->field_vector_local)[ir_dfem] = local_data[ic++];
        }
  }

  VecDestroy(&x);
  VecDestroy(&x_count);

}

void mcpartra::ResidualSourceA::MakeFFQ0Discontinuous()
{
  typedef SpatialDiscretization_PWLD PWLD;
  //======================================== Check correct ff type
  if (resid_ff->spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::ResidualSource3: "
      << "Only PWLD spatial discretizations are supported for now.";
    exit(EXIT_FAILURE);
  }

  const auto& pwl =
    std::dynamic_pointer_cast<PWLD>(resid_ff->spatial_discretization);
  auto& uk_man = resid_ff->unknown_manager;
  auto& field = (*resid_ff->field_vector_local);

  for (const auto& cell : grid->local_cells)
  {
    const size_t num_nodes = pwl->GetCellNumNodes(cell);

    for (size_t m=0; m<uk_man.unknowns.size(); ++m)
      for (size_t g=0; g<uk_man.unknowns[m].num_components; ++g)
      {
        double nodal_avg = 0.0;
        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t dof_map = pwl->MapDOFLocal(cell, i, uk_man, m, g);

          nodal_avg += field[dof_map];
        }
        nodal_avg /= num_nodes;

        for (size_t i=0; i<num_nodes; ++i)
        {
          const int64_t dof_map = pwl->MapDOFLocal(cell, i, uk_man, m, g);

          field[dof_map] = nodal_avg;
        }
      }
  }//for cell

}