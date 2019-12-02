//module:Monte-carlo Utilities
RegisterFunction(chiMonteCarlonCreateSolver);
RegisterFunction(chiMonteCarlonCreateSource);
RegisterFunction(chiMonteCarlonInitialize);
RegisterFunction(chiMonteCarlonExecute);

RegisterNamespace(MCSrcTypes);
AddNamedConstantToNamespace(VACUUM,                1,MCSrcTypes);
AddNamedConstantToNamespace(POINT_SRC,             1,MCSrcTypes);
AddNamedConstantToNamespace(BNDRY_SRC,             2,MCSrcTypes);
AddNamedConstantToNamespace(LOGICVOL_SRC,          3,MCSrcTypes);
AddNamedConstantToNamespace(RESID_SRC,             4,MCSrcTypes);
AddNamedConstantToNamespace(RESID_SRC_SU,          5,MCSrcTypes);
AddNamedConstantToNamespace(RESID_MOC,             6,MCSrcTypes);
AddNamedConstantToNamespace(RESID_MOC_SU,          7,MCSrcTypes);
AddNamedConstantToNamespace(RESIDUAL,              8,MCSrcTypes);
AddNamedConstantToNamespace(RESIDUAL_SU,           9,MCSrcTypes);

RegisterFunction(chiMonteCarlonSetProperty);
RegisterNamespace(MCProperties);
  AddNamedConstantToNamespace(NUM_PARTICLES,               1,MCProperties);
  AddNamedConstantToNamespace(TFC_UPDATE_INTVL,            2,MCProperties);
  AddNamedConstantToNamespace(MONOENERGETIC,               3,MCProperties);
  AddNamedConstantToNamespace(SCATTERING_ORDER,            4,MCProperties);
  AddNamedConstantToNamespace(FORCE_ISOTROPIC,             5,MCProperties);
  AddNamedConstantToNamespace(GROUP_BOUNDS,                6,MCProperties);
  AddNamedConstantToNamespace(TALLY_MERGE_INTVL,           7,MCProperties);
  AddNamedConstantToNamespace(TALLY_MULTIPLICATION_FACTOR, 8,MCProperties);
  AddNamedConstantToNamespace(MAKE_PWLD_SOLUTION,          9,MCProperties);
  AddNamedConstantToNamespace(UNCOLLIDED_ONLY,            10,MCProperties);