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
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,4, 1)
pquad2 = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,5, 5)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)

--############################################### Set boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi;
-- chiLBSSetProperty(phys1,BOUNDARY_CONDITION,ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)

--############################################### Initialize and Execute SteadyStateSolver
chiLBSInitialize(phys1)
single_file = true
is_flux_file = true
chiLBSReadFluxMoments(phys1,"ZMoments",single_file,is_flux_file)
chiLBSSetProperty(phys1,USE_SOURCE_MOMENTS,true)
chiLBSExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)

chiExportFieldFunctionToVTKG(fflist[1],"ZPhiLBS2")
