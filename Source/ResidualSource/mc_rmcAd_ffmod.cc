#include "mc_rmcA_source.h"

#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwl.h"
#include "ChiMath/SpatialDiscretization/FiniteElement/PiecewiseLinear/pwlc.h"
#include "ChiMath/PETScUtils/petsc_utils.h"

#include "chi_log.h"
extern ChiLog& chi_log;

void mcpartra::ResidualSourceA::RemoveFFDiscontinuities()
{
//  size_t num_verts = resid_ff->grid->vertices.size();
//
//  //======================================== Contribute averages
//  std::vector<double> vertex_totals(num_verts,0.0);
//  std::vector<int>    vertex_count(num_verts,0);
//  for (auto& cell : resid_ff->grid->local_cells)
//  {
//    int rmap = (*resid_ff->local_cell_dof_array_address)[cell.local_id];
//    for (int dof=0; dof<cell.vertex_ids.size(); ++dof)
//    {
//      int egrp = 0;
//      int ir = rmap +
//               dof * resid_ff->num_components * resid_ff->num_sets +
//               resid_ff->num_components * 0 +
//               egrp;
//      double phi = (*resid_ff->field_vector_local)[ir];
//
//      vertex_totals[cell.vertex_ids[dof]] += phi;
//      vertex_count[cell.vertex_ids[dof]] += 1;
//    }//for v
//  }//for c
//
//  //======================================== Compute average
//  for (int v=0; v<num_verts; ++v)
//    vertex_totals[v] /= std::max(1,vertex_count[v]);
//
//  //======================================== Reassign field function
//  for (auto& cell : resid_ff->grid->local_cells)
//  {
//    int rmap = (*resid_ff->local_cell_dof_array_address)[cell.local_id];
//    for (int dof=0; dof<cell.vertex_ids.size(); ++dof)
//    {
//      int egrp = 0;
//      int ir = rmap +
//               dof * resid_ff->num_components * resid_ff->num_sets +
//               resid_ff->num_components * 0 +
//               egrp;
//      double phi = vertex_totals[cell.vertex_ids[dof]];
//
//      (*resid_ff->field_vector_local)[ir] = phi;
//    }//for v
//  }//for c




  //======================================== Check correct ff type
  if (resid_ff->spatial_discretization->type !=
      chi_math::SpatialDiscretizationType::PIECEWISE_LINEAR_DISCONTINUOUS)
  {
    chi_log.Log(LOG_ALLERROR)
      << "chi_montecarlon::ResidualSource3: "
      << "Only PWLD spatial discretizations are supported for now.";
    exit(EXIT_FAILURE);
  }

  auto pwl    =
    std::dynamic_pointer_cast<SpatialDiscretization_PWLD>(resid_ff->spatial_discretization);
  auto uk_man = resid_ff->unknown_manager;

  auto pwl_cfem = SpatialDiscretization_PWLC::New(grid,
                    chi_math::finite_element::COMPUTE_UNIT_INTEGRALS);

  auto cfem_num_local_dofs = pwl_cfem->GetNumLocalDOFs(uk_man);
  auto cfem_num_globl_dofs = pwl_cfem->GetNumGlobalDOFs(uk_man);

  Vec x = chi_math::PETScUtils::CreateVector(static_cast<int64_t>(cfem_num_local_dofs),
                                             static_cast<int64_t>(cfem_num_globl_dofs));
  VecSet(x,0.0);
  Vec x_count;
  VecDuplicate(x,&x_count);

  std::vector<int64_t> global_cfem_dofs;
  for (auto& cell : grid->local_cells)
  {
    int i=-1;
    for (uint64_t vid : cell.vertex_ids)
    {
      ++i;
      int64_t ir_dfem = pwl->MapDOFLocal(cell,i,uk_man,0,0);
      int64_t ir_cfem = pwl_cfem->MapDOF(cell,i);

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

  int ic=-1;
  for (auto& cell : grid->local_cells)
  {
    int i=-1;
    for (uint64_t vid : cell.vertex_ids)
    {
      ++i;
      ++ic;
      int64_t ir_dfem = pwl->MapDOFLocal(cell,i,uk_man,0,0);

      (*resid_ff->field_vector_local)[ir_dfem] = local_data[ic];
    }
  }

  VecDestroy(&x);
  VecDestroy(&x_count);

}