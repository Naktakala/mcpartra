#include "sdsolver.h"

#include "chi_log.h"
#include "chi_mpi.h"


//###################################################################
/**Writes a run tape to file.*/
void mcpartra::SourceDrivenSolver::WriteRunTape(const std::string &file_base_name)
{
  auto& chi_log = ChiLog::GetInstance();
  auto& chi_mpi = ChiMPI::GetInstance();

  std::string file_name =
    file_base_name + std::to_string(chi_mpi.location_id) + ".r";

  //============================================= Open file
  chi_log.Log() << "Writing RunTape file " << file_name;
  std::ofstream file(file_name,
                     std::ofstream::binary | //binary file
                     std::ofstream::out |    //no accidental reading
                     std::ofstream::trunc);  //clear file contents when opened

  //============================================= Check file is open
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLWARNING)
      << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Write header
  std::string header_info =
    "Chi-Tech LinearBoltzmann: Flux moments file\n"
    "Header size: 320 bytes\n"
    "Structure(type-info):\n"
    "size_t-num_local_nodes\n"
    "size_t-num_moments\n"
    "size_t-num_groups\n"
    "size_t-num_records\n"
    "Each record:\n"
    "size_t-cell_global_id\n"
    "unsigned int-node_number\n"
    "unsigned int-moment_num\n"
    "unsigned int-group_num\n"
    "double-flux_moment_value\n";

  int header_size = (int)header_info.length();

  char header_bytes[320];
  memset(header_bytes, '-', 320);
  strncpy(header_bytes, header_info.c_str(),std::min(header_size,319));
  header_bytes[319]='\0';

  file << header_bytes;

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;

  size_t num_fv_local_nodes = fv->GetNumLocalDOFs(NODES_ONLY);
  size_t num_fv_local_dofs  = fv->GetNumLocalDOFs(uk_man_fv);

  size_t num_pwl_local_nodes = pwl->GetNumLocalDOFs(NODES_ONLY);
  size_t num_pwl_local_dofs  = pwl->GetNumLocalDOFs(uk_man_pwld);

  size_t num_fv_tallies  = fv_tallies.size();
  size_t num_pwl_tallies = pwl_tallies.size();

  //============================================= Write num_ quantities
  file.write((char*)&num_moments        ,sizeof(size_t));
  file.write((char*)&num_groups         ,sizeof(size_t));
  file.write((char*)&num_fv_local_nodes ,sizeof(size_t));
  file.write((char*)&num_fv_local_dofs  ,sizeof(size_t));
  file.write((char*)&num_pwl_local_nodes,sizeof(size_t));
  file.write((char*)&num_pwl_local_dofs ,sizeof(size_t));
  file.write((char*)&num_fv_tallies     ,sizeof(size_t));
  file.write((char*)&num_pwl_tallies    ,sizeof(size_t));
  file.write((char*)&nps_global         ,sizeof(size_t));

  //============================================= Write per dof data
  try
  {
    for (unsigned int t : fv_tallies)
    {
      size_t dof_counter = 0;
      const auto& tally = grid_tally_blocks[t];
      for (const auto& cell : grid->local_cells)
        for (unsigned int i=0; i<fv->GetCellNumNodes(cell); ++i)
          for (unsigned int m=0; m<num_moments; ++m)
            for (unsigned int g=0; g<num_groups; ++g)
            {
              uint64_t dof_map = fv->MapDOFLocal(cell,i,uk_man_fv,m,g);
              double value     = tally.tally_global[dof_map];
              double value_sqr = tally.tally_sqr_global[dof_map];

              file.write((char*)&cell.global_id,sizeof(uint64_t));
              file.write((char*)&i             ,sizeof(unsigned int));
              file.write((char*)&m             ,sizeof(unsigned int));
              file.write((char*)&g             ,sizeof(unsigned int));
              file.write((char*)&value         ,sizeof(double));
              file.write((char*)&value_sqr     ,sizeof(double));
              ++dof_counter;
            }
      if (dof_counter != num_fv_local_dofs)
        chi_log.Log() << "Dof-counter integrity test failed for FV tally.";
    }//for tally t
    for (unsigned int t : pwl_tallies)
    {
      size_t dof_counter = 0;
      const auto& tally = grid_tally_blocks[t];
      for (const auto& cell : grid->local_cells)
        for (unsigned int i=0; i<pwl->GetCellNumNodes(cell); ++i)
          for (unsigned int m=0; m<num_moments; ++m)
            for (unsigned int g=0; g<num_groups; ++g)
            {
              uint64_t dof_map = pwl->MapDOFLocal(cell,i,uk_man_pwld,m,g);
              double value     = tally.tally_global[dof_map];
              double value_sqr = tally.tally_sqr_global[dof_map];

              file.write((char*)&cell.global_id,sizeof(uint64_t));
              file.write((char*)&i             ,sizeof(unsigned int));
              file.write((char*)&m             ,sizeof(unsigned int));
              file.write((char*)&g             ,sizeof(unsigned int));
              file.write((char*)&value         ,sizeof(double));
              file.write((char*)&value_sqr     ,sizeof(double));
              ++dof_counter;
            }
      if (dof_counter != num_pwl_local_dofs)
        chi_log.Log() << "Dof-counter integrity test failed for PWL tally.";
    }//for tally t
  }
  catch (const std::out_of_range& e)
  {
    chi_log.Log(LOG_ALLWARNING) << __FUNCTION__ << ": array access error.";
  }

  //============================================= Clean-up
  file.close();
}


//###################################################################
/**Reads a run-tape file.*/
void mcpartra::SourceDrivenSolver::ReadRunTape(const std::string &file_name)
{
  auto& chi_log = ChiLog::GetInstance();
  auto& chi_mpi = ChiMPI::GetInstance();

  //============================================= Open file
  chi_log.Log() << "Reading RunTape file " << file_name;
  std::ifstream file(file_name,
                     std::ofstream::binary | //binary file
                     std::ofstream::in);     //no accidental writing

  //============================================= Check file is open
  if (not file.is_open())
  {
    chi_log.Log(LOG_ALLWARNING)
      << __FUNCTION__ << "Failed to open " << file_name;
    return;
  }

  //============================================= Read header
  char header_bytes[320]; header_bytes[319] = '\0';
  file.read(header_bytes,319);

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;

  size_t num_fv_local_nodes = fv->GetNumLocalDOFs(NODES_ONLY);
  size_t num_fv_local_dofs  = fv->GetNumLocalDOFs(uk_man_fv);

  size_t num_pwl_local_nodes = pwl->GetNumLocalDOFs(NODES_ONLY);
  size_t num_pwl_local_dofs  = pwl->GetNumLocalDOFs(uk_man_pwld);

  size_t num_fv_tallies  = fv_tallies.size();
  size_t num_pwl_tallies = pwl_tallies.size();

  //============================================= Define file quantities
  size_t file_num_moments         = 0;
  size_t file_num_groups          = 0;
  size_t file_num_fv_local_nodes  = 0;
  size_t file_num_fv_local_dofs   = 0;
  size_t file_num_pwl_local_nodes = 0;
  size_t file_num_pwl_local_dofs  = 0;
  size_t file_num_fv_tallies      = 0;
  size_t file_num_pwl_tallies     = 0;
  size_t file_nps_global          = 0;

  //============================================= Read file quantities
  file.read((char*)&file_num_moments         , sizeof(size_t));
  file.read((char*)&file_num_groups          , sizeof(size_t));
  file.read((char*)&file_num_fv_local_nodes  , sizeof(size_t));
  file.read((char*)&file_num_fv_local_dofs   , sizeof(size_t));
  file.read((char*)&file_num_pwl_local_nodes , sizeof(size_t));
  file.read((char*)&file_num_pwl_local_dofs  , sizeof(size_t));
  file.read((char*)&file_num_fv_tallies      , sizeof(size_t));
  file.read((char*)&file_num_pwl_tallies     , sizeof(size_t));
  file.read((char*)&file_nps_global          , sizeof(size_t));

  //============================================= Check compatibility
  if (file_num_moments         != num_moments         or
      file_num_groups          != num_groups          or
      file_num_fv_local_nodes  != num_fv_local_nodes  or
      file_num_fv_local_dofs   != num_fv_local_dofs   or
      file_num_pwl_local_nodes != num_pwl_local_nodes or
      file_num_pwl_local_dofs  != num_pwl_local_dofs  or
      file_num_fv_tallies      != num_fv_tallies      or
      file_num_pwl_tallies     != num_pwl_tallies)
  {
    std::stringstream outstr;

    outstr <<
    "file_num_moments         " << file_num_moments         << " vs " << file_num_moments         << "\n"
    "file_num_groups          " << file_num_groups          << " vs " << file_num_groups          << "\n"
    "file_num_fv_local_nodes  " << file_num_fv_local_nodes  << " vs " << file_num_fv_local_nodes  << "\n"
    "file_num_fv_local_dofs   " << file_num_fv_local_dofs   << " vs " << file_num_fv_local_dofs   << "\n"
    "file_num_pwl_local_nodes " << file_num_pwl_local_nodes << " vs " << file_num_pwl_local_nodes << "\n"
    "file_num_pwl_local_dofs  " << file_num_pwl_local_dofs  << " vs " << file_num_pwl_local_dofs  << "\n"
    "file_num_fv_tallies      " << file_num_fv_tallies      << " vs " << file_num_fv_tallies      << "\n"
    "file_num_pwl_tallies     " << file_num_pwl_tallies     << " vs " << file_num_pwl_tallies     << "\n";

    chi_log.Log(LOG_ALLWARNING)
      << "Incompatible data found in runtape file " << file_name << "\n"
      << outstr.str();
    file.close();
    return;
  }//if not compatible

  //============================================= Commit to reading the file
  nps_global += file_nps_global;
  try
  {
    for (unsigned int t : fv_tallies)
    {
      size_t dof_counter = 0;
      auto& tally = grid_tally_blocks[t];
      for (size_t dof=0; dof < num_fv_local_dofs; ++dof)
      {
        uint64_t     cell_global_id = 0;
        unsigned int node = 0;
        unsigned int moment = 0;
        unsigned int group = 0;
        double       value = 0.0;
        double       value_sqr = 0.0;

        file.read((char*)&cell_global_id,sizeof(uint64_t));
        file.read((char*)&node          ,sizeof(unsigned int));
        file.read((char*)&moment        ,sizeof(unsigned int));
        file.read((char*)&group         ,sizeof(unsigned int));
        file.read((char*)&value         ,sizeof(double));
        file.read((char*)&value_sqr     ,sizeof(double));

        const auto& cell = grid->cells[cell_global_id];

        size_t dof_map = fv->MapDOFLocal(cell,node,uk_man_fv,moment,group);

        tally.tally_global    [dof_map] += value;
        tally.tally_sqr_global[dof_map] += value_sqr;
        ++dof_counter;
      }//for dof
      if (dof_counter != num_fv_local_dofs)
        chi_log.Log() << "Dof-counter integrity test failed for FV tally.";
    }//for tally t
    for (unsigned int t : pwl_tallies)
    {
      size_t dof_counter = 0;
      auto& tally = grid_tally_blocks[t];
      for (size_t dof=0; dof < num_pwl_local_dofs; ++dof)
      {
        uint64_t     cell_global_id = 0;
        unsigned int node = 0;
        unsigned int moment = 0;
        unsigned int group = 0;
        double       value = 0.0;
        double       value_sqr = 0.0;

        file.read((char*)&cell_global_id,sizeof(uint64_t));
        file.read((char*)&node          ,sizeof(unsigned int));
        file.read((char*)&moment        ,sizeof(unsigned int));
        file.read((char*)&group         ,sizeof(unsigned int));
        file.read((char*)&value         ,sizeof(double));
        file.read((char*)&value_sqr     ,sizeof(double));

        const auto& cell = grid->cells[cell_global_id];

        size_t dof_map = pwl->MapDOFLocal(cell,node,uk_man_pwld,moment,group);

        tally.tally_global    [dof_map] += value;
        tally.tally_sqr_global[dof_map] += value_sqr;
        ++dof_counter;
      }//for dof
      if (dof_counter != num_pwl_local_dofs)
        chi_log.Log() << "Dof-counter integrity test failed for PWL tally.";
    }//for tally t
  }
  catch (const std::out_of_range& e)
  {
    chi_log.Log(LOG_ALLWARNING) << __FUNCTION__ << ": array access error.";
  }


  file.close();
}

//###################################################################
/**Write LBS-formatted source moments.*/
void mcpartra::SourceDrivenSolver::WriteLBSFluxMoments(const std::string &file_name)
{
  auto& chi_log = ChiLog::GetInstance();
  auto& chi_mpi = ChiMPI::GetInstance();

  //============================================= Get relevant items
  auto NODES_ONLY = ChiMath::UNITARY_UNKNOWN_MANAGER;
  uint64_t num_local_nodes = pwl->GetNumGlobalDOFs(NODES_ONLY);
  uint64_t num_local_dofs  = pwl->GetNumGlobalDOFs(uk_man_pwld);
  uint64_t num_local_cells = grid->GetGlobalNumberOfCells();

  //============================================= Open file sequentially on
  //                                              each processor
  chi_log.Log() << "Writing flux-moments to " << file_name << ".";
  std::ofstream file;
  auto mode_all = std::ofstream::binary | std::ofstream::out;

  //================================================== Write size quantities
  //                                                   and cell nodal info
  for (int location=0; location < chi_mpi.process_count; ++location)
  {
    MPI_Barrier(MPI_COMM_WORLD);

    //=========================================== If location does not have scope
    //                                            it will go wait at the barrier
    if (location != chi_mpi.location_id) continue;

    //=========================================== Home location wipes file,
    //                                            others append to it
    if (location == 0 and chi_mpi.location_id == location)
      file.open(file_name, mode_all | std::ofstream::trunc);
    else
      file.open(file_name, mode_all | std::ofstream::app);

    //=========================================== Check file is open
    if (not file.is_open())
    {
      chi_log.Log(LOG_ALLWARNING)
        << __FUNCTION__ << "Failed to open " << file_name;
      continue;
    }

    //=========================================== Home location writes header
    if (location == 0 and chi_mpi.location_id == location)
    {
      std::string header_info =
        "Chi-Tech LinearBoltzmann: Flux moments file\n"
        "Header size: 500 bytes\n"
        "Structure(type-info):\n"
        "uint64_t num_local_nodes\n"
        "uint64_t num_moments\n"
        "uint64_t num_groups\n"
        "uint64_t num_records\n"
        "uint64_t num_cells\n"
        "Each cell:\n"
        "  uint64_t cell_global_id\n"
        "  uint64_t   num_nodes\n"
        "  Each node:\n"
        "    double   x_position\n"
        "    double   y_position\n"
        "    double   z_position\n"
        "Each record:\n"
        "  uint64_t     cell_global_id\n"
        "  unsigned int node_number\n"
        "  unsigned int moment_num\n"
        "  unsigned int group_num\n"
        "  double       flux_moment_value\n";

      int header_size = (int) header_info.length();

      char header_bytes[500];
      memset(header_bytes, '-', 500);
      strncpy(header_bytes, header_info.c_str(), std::min(header_size, 499));
      header_bytes[499] = '\0';

      file << header_bytes;

      //============================================= Write num_ quantities
      uint64_t num_moments_t = num_moments;
      uint64_t num_groups_t  = num_groups;
      file.write((char*)&num_local_nodes,sizeof(uint64_t));
      file.write((char*)&num_moments_t  ,sizeof(uint64_t));
      file.write((char*)&num_groups_t   ,sizeof(uint64_t));
      file.write((char*)&num_local_dofs ,sizeof(uint64_t));
      file.write((char*)&num_local_cells,sizeof(uint64_t));
    }

    //============================================= Write nodal positions for
    //                                              each cell
    for (const auto& cell : grid->local_cells)
    {
      file.write((char *) &cell.global_id, sizeof(uint64_t));

      uint64_t num_nodes = pwl->GetCellNumNodes(cell);
      file.write((char *) &num_nodes, sizeof(uint64_t));

      auto   node_locations = pwl->GetCellNodeLocations(cell);
      for (const auto& node : node_locations)
      {
        file.write((char *) &node.x, sizeof(double));
        file.write((char *) &node.y, sizeof(double));
        file.write((char *) &node.z, sizeof(double));
      }//for node
    }//for cell

    file.close();
  }//for each location

  //================================================== Write dof data
  for (int location=0; location < chi_mpi.process_count; ++location)
  {
    MPI_Barrier(MPI_COMM_WORLD);

    //=========================================== If location does not have scope
    //                                            it will go wait at the barrier
    if (location != chi_mpi.location_id) continue;

    //=========================================== All locations just append
    file.open(file_name, mode_all | std::ofstream::app);

    //=========================================== Check file is open
    if (not file.is_open())
    {
      chi_log.Log(LOG_ALLWARNING)
        << __FUNCTION__ << "Failed to open " << file_name << " again.";
      continue;
    }

    //============================================= Write per dof data
    auto& pwl_tally = grid_tally_blocks.at(pwl_tallies.front()).tally_global;
    try {
      for (const auto& cell : grid->local_cells)
        for (unsigned int i=0; i<pwl->GetCellNumNodes(cell); ++i)
          for (unsigned int m=0; m<num_moments; ++m)
            for (unsigned int g=0; g<num_groups; ++g)
            {
              uint64_t dof_map = pwl->MapDOFLocal(cell,i,uk_man_pwld,m,g);
              double value = pwl_tally.at(dof_map);

              file.write((char*)&cell.global_id,sizeof(uint64_t));
              file.write((char*)&i             ,sizeof(unsigned int));
              file.write((char*)&m             ,sizeof(unsigned int));
              file.write((char*)&g             ,sizeof(unsigned int));
              file.write((char*)&value         ,sizeof(double));
            }
    }
    catch (const std::out_of_range& e)
    {
      chi_log.Log(LOG_ALLWARNING) << __FUNCTION__ << ": The given flux_moments-"
                                  << "vector was accessed out of range.";
    }


    file.close();
  }//for each location

}

//###################################################################
/**Reads a binary cell-wise importance map associated with a certain QOI.*/
void mcpartra::SourceDrivenSolver::ReadImportanceMap(const std::string &file_name)
{
  const std::string fname = __FUNCTION__;
  auto& chi_log = ChiLog::GetInstance();

  chi_log.Log() << "MCParTra: Reading importance map.";

  //============================================= Check grid is available
  if (grid == nullptr)
    throw std::logic_error(fname + ": Grid not available yet. This function"
                                   " must be called after initialization.");

  //============================================= Open file
  std::ifstream file(file_name, std::ios_base::in | std::ios_base::binary);

  //============================================= Check file is open
  if (not file.is_open())
    throw std::logic_error(fname + ": Could not open file " + file_name);

  //============================================= Read header
  char header_bytes[400]; header_bytes[399] = '\0';
  file.read(header_bytes,399);

  uint64_t file_num_global_cells;
  uint64_t file_num_groups      ;
  uint64_t file_num_records     ;

  file.read((char*)&file_num_global_cells, sizeof(uint64_t));
  file.read((char*)&file_num_groups      , sizeof(uint64_t));
  file.read((char*)&file_num_records     , sizeof(uint64_t));

  //============================================= Error check
  size_t sim_num_global_cells = grid->GetGlobalNumberOfCells();
  if (file_num_global_cells != sim_num_global_cells)
    throw std::logic_error(fname + ": Data contained in " + file_name +
    " is not compatible with the current simulation.  The cell count is "
    "not equal " + std::to_string(file_num_global_cells) + " (file) vs " +
    std::to_string(sim_num_global_cells) + " (sim).");

  if (file_num_groups != num_groups)
    throw std::logic_error(fname + ": Data contained in " + file_name +
    " is not compatible with the current simulation. The number of groups are "
    "not equal " + std::to_string(file_num_groups) + " (file) vs " +
    std::to_string(num_groups) + " (sim).");


  //============================================= Size importance data structs
  size_t num_local_cells = grid->local_cells.size();
  importance_num_groups = file_num_groups;
  typedef chi_mesh::Vector3 vec3;

  size_t num_importances = num_local_cells * file_num_groups;
  local_cell_importance.assign(num_importances, 1.0);
  local_cell_importance_setting.assign(num_importances, 1.0);
  local_cell_importance_directions.assign(num_importances, vec3(0.0,0.0,0.0));
  local_cell_importance_exp_coeffs.assign(num_importances, {0.0,0.0});

  //============================================= Read records
  for (uint64_t r=0; r < file_num_records; ++r)
  {
    uint64_t cell_global_id;
    unsigned int group_id;
    std::array<double,4> p1_moments = {0.0,0.0,0.0,0.0};
    double coeff_a;
    double coeff_b;

    file.read((char*)&cell_global_id       , sizeof(uint64_t));
    file.read((char*)&group_id             , sizeof(unsigned int));
    file.read((char*)p1_moments.data()     , sizeof(double)*4);
    file.read((char*)&coeff_a              , sizeof(double));
    file.read((char*)&coeff_b              , sizeof(double));

    if (grid->IsCellLocal(cell_global_id))
    {
      const auto& cell = grid->cells[cell_global_id];

      size_t dof_map = cell.local_id * importance_num_groups + group_id;

      double phi = p1_moments[0];
      chi_mesh::Vector3 omega_J(p1_moments[1],
                                p1_moments[2],
                                p1_moments[3]);
      omega_J.Normalize();

      local_cell_importance[dof_map]            = phi;
      local_cell_importance_setting[dof_map]    = phi;
      local_cell_importance_directions[dof_map] = omega_J;
      local_cell_importance_exp_coeffs[dof_map] = std::make_pair(coeff_a,coeff_b);
    }
  }//for record

  file.close();

  chi_log.Log() << "MCParTra: Done reading importance map.";

  for (auto& src : sources)
    src->BiasCDFs(false);

  //============================================= Make unknown managers
  chi_math::UnknownManager uk_man;
  uk_man.AddUnknown(chi_math::UnknownType::VECTOR_N,3);

  uk_man.SetUnknownTextName(0, "Comp");

  //============================================= Make vectors compatible with
  //                                              field functions
  std::vector<double> data_omega_J;
  std::vector<double> data_exp_coeffs;

  {
    auto& data_vector = data_omega_J;
    size_t num_local_dofs = fv->GetNumLocalDOFs(uk_man);

    data_vector.assign(num_local_dofs, 0.0);

    uint64_t cell_local_id = 0;
    for (const auto& omega_J : local_cell_importance_directions)
    {
      const auto& cell = grid->local_cells[cell_local_id];
      int64_t dof_map = fv->MapDOFLocal(cell, 0, uk_man, 0, 0);

      data_vector[dof_map+0] = omega_J.x;
      data_vector[dof_map+1] = omega_J.y;
      data_vector[dof_map+2] = omega_J.z;

      ++cell_local_id;
    }
  }
  {
    auto& data_vector = data_exp_coeffs;
    size_t num_local_dofs = fv->GetNumLocalDOFs(uk_man);

    data_vector.assign(num_local_dofs, 0.0);

    uint64_t cell_local_id = 0;
    for (const auto& coeffs : local_cell_importance_exp_coeffs)
    {
      const auto& cell = grid->local_cells[cell_local_id];
      int64_t dof_map = fv->MapDOFLocal(cell, 0, uk_man, 0, 0);

      data_vector[dof_map+0] = coeffs.first;
      data_vector[dof_map+1] = coeffs.second;

      ++cell_local_id;
    }
  }

  //============================================= Export interior source
  //                                              as FieldFunction
  auto fv_sd = std::dynamic_pointer_cast<SpatialDiscretization>(fv);
  auto R_ff = std::make_shared<chi_physics::FieldFunction>(
    "omega_J",                                 //Text name
    fv_sd,                                     //Spatial Discretization
    &data_omega_J,                             //Data
    uk_man,                                  //Nodal variable structure
    0, 0);                                     //Reference variable and component

  R_ff->ExportMultipleFFToVTK("Z_J", {R_ff});
}