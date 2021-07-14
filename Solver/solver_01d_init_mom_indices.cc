#include "solver_montecarlon.h"

#include "chi_log.h"
extern ChiLog& chi_log;

#include "chi_mpi.h"
extern ChiMPI& chi_mpi;

//###################################################################
/**Initialize moment indices.*/
void mcpartra::Solver::InitMomentIndices()
{
  const auto& first_cell = grid->local_cells[0];

  m_to_ell_em_map.clear();
  if (first_cell.Type() == chi_mesh::CellType::SLAB)
    for (int ell=0; ell<=options.scattering_order; ell++)
      m_to_ell_em_map.emplace_back(ell,0);
  else if (first_cell.Type() == chi_mesh::CellType::POLYGON)
    for (int ell=0; ell<=options.scattering_order; ell++)
      for (int m=-ell; m<=ell; m+=2)
      {
        if (ell == 0 or m != 0)
          m_to_ell_em_map.emplace_back(ell,m);
      }
  else if (first_cell.Type() == chi_mesh::CellType::POLYHEDRON)
    for (int ell=0; ell<=options.scattering_order; ell++)
      for (int m=-ell; m<=ell; m++)
        m_to_ell_em_map.emplace_back(ell,m);

  num_moments = m_to_ell_em_map.size();
}

