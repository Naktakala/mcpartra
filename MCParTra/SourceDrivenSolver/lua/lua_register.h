//module:Monte-carlo Utilities
RegisterFunction(chiMonteCarlonCreateSolver);
RegisterFunction(chiMonteCarlonCreateSource);
RegisterFunction(chiMonteCarlonInitialize);
RegisterFunction(chiMonteCarlonExecute);

RegisterFunction(chiMonteCarlonReadRuntape);
RegisterFunction(chiMonteCarlonWriteLBSFluxMoments);

RegisterNamespace(MCSrcTypes);
AddNamedConstantToNamespace(POINT_SRC,             1,MCSrcTypes);
AddNamedConstantToNamespace(BNDRY_SRC,             2,MCSrcTypes);
AddNamedConstantToNamespace(MATERIAL_SRC,          3,MCSrcTypes);
AddNamedConstantToNamespace(RESIDUAL_TYPE_A,       4,MCSrcTypes);
AddNamedConstantToNamespace(RESIDUAL_TYPE_B,       5,MCSrcTypes);

RegisterFunction(chiMonteCarlonSetImportances)
RegisterFunction(chiMonteCarlonSetProperty2);
RegisterNamespace(MCProperties);

RegisterFunction(chiMonteCarlonAddCustomVolumeTally);