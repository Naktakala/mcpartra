-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=0.50758 and 2.52527e-04
num_procs = 1





--############################################### Check num_procs
if (check_num_procs==nil and chi_number_of_processes ~= num_procs) then
    chiLog(LOG_0ERROR,"Incorrect amount of processors. " ..
                      "Expected "..tostring(num_procs)..
                      ". Pass check_num_procs=false to override if possible.")
    os.exit(false)
end

--############################################### Setup mesh
chiMeshHandlerCreate()

chiUnpartitionedMeshFromEnsightGold("RegressionTests/BLOCK_CITYV2a.case")

region1 = chiRegionCreate()
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute();

--############################################### Set Material IDs
-- vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
-- chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

chiRegionExportMeshToVTK(region1, "ZMesh")

--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");
materials[3] = chiPhysicsAddMaterial("Test Material3");
materials[4] = chiPhysicsAddMaterial("Test Material4");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[3],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[4],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[3],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[4],ISOTROPIC_MG_SOURCE)

num_groups = 116
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"RegressionTests/xs_air_50rh_116La.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"RegressionTests/xs_ss304_116.data")
chiPhysicsMaterialSetProperty(materials[3],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"RegressionTests/xs_ss304_116.data")
chiPhysicsMaterialSetProperty(materials[4],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"RegressionTests/xs_air_50rh_116La.data")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
-- src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
-- src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
-- src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[3],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[4],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
if (seed == nil) then seed=chi_location_id end
phys1 = chiMonteCarlonCreateSolver(seed)
chiSolverAddRegion(phys1,region1)

-- chiMonteCarlonCreateSource(phys1,MCSrcTypes.BNDRY_SRC,1);
chiMonteCarlonCreateSource(phys1,MCSrcTypes.MATERIAL_SRC);

-- chiMonteCarlonSetProperty(phys1,MCProperties.NUM_PARTICLES,1e6)
-- chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MERGE_INTVL,1e5)
-- chiMonteCarlonSetProperty(phys1,MCProperties.SCATTERING_ORDER,0)
-- chiMonteCarlonSetProperty(phys1,MCProperties.MONOENERGETIC,false)
-- chiMonteCarlonSetProperty(phys1,MCProperties.FORCE_ISOTROPIC,false)
-- chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0)
-- chiMonteCarlonSetProperty(phys1,MCProperties.MAKE_PWLD_SOLUTION,true)
-- chiMonteCarlonSetProperty(phys1,MCProperties.UNCOLLIDED_ONLY,true)

chiMonteCarlonSetProperty2(phys1,"NUM_PARTICLES"              ,0.1e5)
chiMonteCarlonSetProperty2(phys1,"TALLY_MERGE_INTVL"          ,1e4)
chiMonteCarlonSetProperty2(phys1,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys1,"MONOENERGETIC"              ,false)
chiMonteCarlonSetProperty2(phys1,"FORCE_ISOTROPIC"            ,false)
chiMonteCarlonSetProperty2(phys1,"TALLY_MULTIPLICATION_FACTOR",1.0)
chiMonteCarlonSetProperty2(phys1,"MAKE_PWLD_SOLUTION"         ,true)
chiMonteCarlonSetProperty2(phys1,"UNCOLLIDED_ONLY"            ,true)
if (run_tape_basename == nil) then run_tape_basename="ZRunTape" end
chiMonteCarlonSetProperty2(phys1,"RUN_TAPE_BASE_NAME"         ,run_tape_basename)

chiMonteCarlonInitialize(phys1)

chiMonteCarlonReadRuntape(phys1, "ZRunTape00.r")
chiMonteCarlonReadRuntape(phys1, "ZRunTape10.r")
chiMonteCarlonReadRuntape(phys1, "ZRunTape20.r")
chiMonteCarlonReadRuntape(phys1, "ZRunTape30.r")

chiMonteCarlonExecute(phys1)

chiMonteCarlonWriteLBSFluxMoments(phys1, "ZMoments.data");
--
-- --############################################### Setup LBS Physics
-- phys1 = chiLBSCreateSolver()
-- chiSolverAddRegion(phys1,region1)
--
-- --========== Groups
-- grp = {}
-- for g=1,num_groups do
--     grp[g] = chiLBSCreateGroup(phys1)
-- end
--
-- --========== ProdQuad
-- pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,6, 6)
--
-- --========== Groupset def
-- gs0 = chiLBSCreateGroupset(phys1)
-- cur_gs = gs0
-- chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
-- chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
-- chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
-- chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
-- chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
-- chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
-- chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
-- chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
-- chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
-- --chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
-- --chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--
-- --========== Boundary conditions
-- bsrc={}
-- for g=1,num_groups do
--     bsrc[g] = 0.0
-- end
-- bsrc[1] = 1.0/4.0/math.pi
-- -- chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
-- --         YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
-- --chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
-- --        XMIN,LBSBoundaryTypes.REFLECTING);
-- --chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
-- --        XMAX,LBSBoundaryTypes.REFLECTING);
--
-- --========== Solvers
-- chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
-- chiLBSSetProperty(phys1,SCATTERING_ORDER,0)
--
-- chiLBSInitialize(phys1)
-- -- chiLBSExecute(phys1)
--
-- fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

chiExportFieldFunctionToVTKG(0,"ZPhiMC")
chiExportFieldFunctionToVTKG(num_groups,"ZPhiMCPWL")
-- chiExportFieldFunctionToVTKG(fflist[1],"ZPhiLBS")
