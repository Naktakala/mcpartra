#include "mc_material_source.h"

#include "SourceDrivenSolver/sdsolver.h"

#include "ChiPhysics/PhysicsMaterial/material_property_isotropic_mg_src.h"

#include "ChiPhysics/chi_physics.h"
extern ChiPhysics&  chi_physics_handler;

#include "ChiLog/chi_log.h"
extern ChiLog& chi_log;

//###################################################################
/**Initialize material source.*/
void mcpartra::MaterialSource::
  Initialize(chi_mesh::MeshContinuumPtr& ref_grid,
             std::shared_ptr<SpatialDiscretization_FV>& ref_fv_sdm,
             size_t ref_num_groups,
             const std::vector<std::pair<int,int>>& ref_m_to_ell_em_map,
             const std::vector<CellGeometryData>& ref_cell_geometry_info)
{
  chi_log.Log(LOG_0) << "Initializing Material Sources";

  grid = ref_grid;
  fv_sdm = ref_fv_sdm;
  num_groups = ref_num_groups;
  m_to_ell_em_map = ref_m_to_ell_em_map;

  size_t num_materials = chi_physics_handler.material_stack.size();

  //============================================= List cells with mat-sources
  cell_elements.clear();
  matid_has_q_flags.assign(num_materials, false);
  matid_q_map.clear();
  for (auto& cell : ref_grid->local_cells)
  {
    auto mat_id = static_cast<int>(cell.material_id);
    auto mat = chi_physics_handler.material_stack[mat_id];
    for (const auto& prop : mat->properties)
      if (prop->Type() == chi_physics::PropertyType::ISOTROPIC_MG_SOURCE)
      {
        cell_elements[cell.local_id] = std::move(GetCellVolumeElements(cell, grid));
        matid_has_q_flags[mat_id]    = true;
        matid_q_map[mat_id]          = std::static_pointer_cast<IsoMGSrc>(prop);
        break;
      }
  }

  //============================================= Determine group-wise
  //                                              source strengths
  //                                              Unnormalized PDF
  double IntV_Q_total_local = 0.0;
  std::vector<double> IntV_Q_g(num_groups, 0.0);
  group_elements.resize(num_groups);
  {
    for (auto& [cell_local_id, vec_of_elements] : cell_elements)
    {
      const auto& cell = grid->local_cells[cell_local_id];
      const auto fv_view = ref_fv_sdm->MapFeView(cell.local_id);

      const auto mat_id = static_cast<int>(cell.material_id);
      const auto mat = chi_physics_handler.material_stack[mat_id];

      if (matid_has_q_flags[mat_id])
      {
        auto src = matid_q_map[mat_id];
        for (size_t g=0; g<src->source_value_g.size(); ++g)
        {
          double Q_g = src->source_value_g[g];

          if (not (std::fabs(Q_g) > 0.0)) continue;

          IntV_Q_g[g] += fv_view->volume*Q_g;

          for (auto& element : vec_of_elements)
            group_elements[g].emplace_back(Q_g * element.Volume(), element);
        }//for g
      }//if has source
    }//for cell

    for (auto val : IntV_Q_g)
      IntV_Q_total_local += val;
  }

  //======================================== Compute normalized pdf used during
  //                                         biasing
  group_element_pdf.resize(num_groups);
  for (size_t g=0; g<num_groups; ++g)
  {
    size_t num_elems = group_elements[g].size();
    group_element_pdf[g].resize(num_elems, 0.0);
    for (size_t elem=0; elem<num_elems; ++elem)
      group_element_pdf[g][elem] = group_elements[g][elem].first / IntV_Q_g[g];
  }

  //============================================= Set source strength for
  //                                              solver
  {
    double IntV_Q_total_globl = 0.0;
    MPI_Allreduce(&IntV_Q_total_local,     // sendbuf
                  &IntV_Q_total_globl,     // recvbuf
                  1, MPI_DOUBLE,           // count + datatype
                  MPI_SUM,                 // operation
                  MPI_COMM_WORLD);         // communicator

    SetSourceRates(IntV_Q_total_local,IntV_Q_total_globl);
  }

  //============================================= Construct groupwise cdf
  group_cdf.clear();
  group_cdf.resize(num_groups, 0.0);
  {
    double running_total = 0.0;
    for (size_t g=0; g<num_groups; ++g)
    {
      running_total += IntV_Q_g[g];
      group_cdf[g] = (IntV_Q_total_local > 0.0)?
        running_total / IntV_Q_total_local : 0.0;
    }//for g
  }

  //============================================= Within each group, construct
  //                                              source element cdf
  group_element_cdf.clear();
  group_element_cdf.resize(num_groups);
  {
    size_t g = 0;
    for (auto& group_source : group_elements)
    {
      group_element_cdf[g].resize(group_source.size(),0.0);
      double elem_running_total = 0.0;
      size_t elem = 0;
      for (auto& src_element : group_source)
      {
        elem_running_total += src_element.first;
        group_element_cdf[g][elem] = (IntV_Q_g[g] > 0.0)?
          elem_running_total/IntV_Q_g[g] : 0.0;
        ++elem;
      }//for elem
      ++g;
    }//for g
  }

  chi_log.Log(LOG_0) << "Done initializing Material Sources";

  BiasCDFs(false); //This just copies cdfs and corrections
}