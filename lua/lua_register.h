//module:Monte-carlo Utilities
RegisterFunction(chiMonteCarlonCreateSolver)
RegisterFunction(chiMonteCarlonCreateSource)
RegisterFunction(chiMonteCarlonInitialize)
RegisterFunction(chiMonteCarlonExecute)
  RegisterConstant(MC_POINT_SRC,             1);
  RegisterConstant(MC_BNDRY_SRC,             2);
  RegisterConstant(MC_ALL_BOUNDARIES,         -1);
  RegisterConstant(MC_LOGICVOL_SRC,          3);
  RegisterConstant(MC_RESID_SRC,             4);
  RegisterConstant(MC_RESID_SRC_SU,          5);
  RegisterConstant(MC_RESID_MOC,             6);
  RegisterConstant(MC_RESID_MOC_SU,          7);
RegisterFunction(chiMonteCarlonSetProperty)
  RegisterConstant(MC_NUM_PARTICLES,               1);
  RegisterConstant(MC_TFC_UPDATE_INTVL,            2);
  RegisterConstant(MC_MONOENERGETIC,               3);
  RegisterConstant(MC_SCATTERING_ORDER,            4);
  RegisterConstant(MC_FORCE_ISOTROPIC,             5);
  RegisterConstant(MC_GROUP_BOUNDS,                6);
  RegisterConstant(MC_TALLY_MERGE_INTVL,           7);
  RegisterConstant(MC_TALLY_MULTIPLICATION_FACTOR, 8);
  RegisterConstant(MC_MAKE_PWLD_SOLUTION,          9);
  RegisterConstant(MC_UNCOLLIDED_ONLY,             10);