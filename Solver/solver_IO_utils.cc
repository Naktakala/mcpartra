#include "solver_montecarlon.h"

#include "chi_log.h"
#include "chi_mpi.h"

//###################################################################
/**Writes a run tape to file.*/
void mcpartra::Solver::WriteRunTape(const std::string &file_base_name)
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
      const auto& tally = grid_tally_blocks[t];
      for (const auto& cell : grid->local_cells)
        for (unsigned int i=0; i<fv->GetCellNumNodes(cell); ++i)
          for (unsigned int m=0; m<num_moments; ++m)
            for (unsigned int g=0; g<num_groups; ++g)
            {
              uint64_t dof_map = fv->MapDOFLocal(cell,i,uk_man_fv,m,g);
              double value     = tally.tally_global[dof_map];
              double value_sqr = tally.tally_sqr_global[dof_map];

              file.write((char*)&cell.global_id,sizeof(size_t));
              file.write((char*)&i             ,sizeof(unsigned int));
              file.write((char*)&m             ,sizeof(unsigned int));
              file.write((char*)&g             ,sizeof(unsigned int));
              file.write((char*)&value         ,sizeof(double));
              file.write((char*)&value_sqr     ,sizeof(double));
            }
    }//for tally t
    for (unsigned int t : pwl_tallies)
    {
      const auto& tally = grid_tally_blocks[t];
      for (const auto& cell : grid->local_cells)
        for (unsigned int i=0; i<pwl->GetCellNumNodes(cell); ++i)
          for (unsigned int m=0; m<num_moments; ++m)
            for (unsigned int g=0; g<num_groups; ++g)
            {
              uint64_t dof_map = pwl->MapDOFLocal(cell,i,uk_man_pwld,m,g);
              double value     = tally.tally_global[dof_map];
              double value_sqr = tally.tally_sqr_global[dof_map];

              file.write((char*)&cell.global_id,sizeof(size_t));
              file.write((char*)&i             ,sizeof(unsigned int));
              file.write((char*)&m             ,sizeof(unsigned int));
              file.write((char*)&g             ,sizeof(unsigned int));
              file.write((char*)&value         ,sizeof(double));
              file.write((char*)&value_sqr     ,sizeof(double));
            }
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
void mcpartra::Solver::ReadRunTape(const std::string &file_name)
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
      auto& tally = grid_tally_blocks[t];
      for (size_t dof=0; dof < num_fv_local_dofs; ++dof)
      {
        size_t       cell_global_id = 0;
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
      }//for dof
    }//for tally t
    for (unsigned int t : pwl_tallies)
    {
      auto& tally = grid_tally_blocks[t];
      for (size_t dof=0; dof < num_fv_local_dofs; ++dof)
      {
        size_t       cell_global_id = 0;
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

        size_t dof_map = pwl->MapDOFLocal(cell,node,uk_man_fv,moment,group);

        tally.tally_global    [dof_map] += value;
        tally.tally_sqr_global[dof_map] += value_sqr;
      }//for dof
    }//for tally t
  }
  catch (const std::out_of_range& e)
  {
    chi_log.Log(LOG_ALLWARNING) << __FUNCTION__ << ": array access error.";
  }


  file.close();
}

