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

mesh={}
N=10
L=5
xmin = -L/2
dx = L/N
for i=1,(N+1) do
    k=i-1
    mesh[i] = xmin + k*dx
end
chiMeshCreateUnpartitioned3DOrthoMesh(mesh,mesh,mesh)
chiVolumeMesherExecute();

--############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 168
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"RegressionTests/xs_3_170.cxs")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        CHI_XSFILE,"RegressionTests/xs_3_170.cxs")

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--############################################### Setup Physics
phys1 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys1,region1)

-- chiMonteCarlonCreateSource(phys1,MCSrcTypes.BNDRY_SRC,1);
chiMonteCarlonCreateSource(phys1,MCSrcTypes.MATERIAL_SRC);

chiMonteCarlonSetProperty(phys1,MCProperties.NUM_PARTICLES,1e6)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MERGE_INTVL,1e5)
chiMonteCarlonSetProperty(phys1,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys1,MCProperties.MONOENERGETIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.FORCE_ISOTROPIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0)
chiMonteCarlonSetProperty(phys1,MCProperties.MAKE_PWLD_SOLUTION,true)
chiMonteCarlonSetProperty(phys1,MCProperties.UNCOLLIDED_ONLY,true)

chiMonteCarlonInitialize(phys1)
-- os.exit()
chiMonteCarlonExecute(phys1)

--############################################### Setup LBS Physics
phys1 = chiLBSCreateSolver()
chiSolverAddRegion(phys1,region1)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys1)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,16, 8)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_SINGLE)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES_CYCLES)
chiLBSGroupsetSetResidualTolerance(phys1,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys1,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys1,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys1,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys1,cur_gs,30,1.0e-4,false," ")

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/4.0/math.pi
-- chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--         YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        XMIN,LBSBoundaryTypes.REFLECTING);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        XMAX,LBSBoundaryTypes.REFLECTING);

--========== Solvers
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBSInitialize(phys1)
chiLBSExecute(phys1)

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)


--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,-0.01,-L/2, 0.0)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,-0.01,-L/2, 0.0)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

-- for k=1,2 do
--     chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,k-1)
-- end
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,num_groups)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[1])
--chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist[2])


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)


--


if (chi_location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end

chiExportFieldFunctionToVTKG(0,"ZPhiMC")
chiExportFieldFunctionToVTKG(num_groups,"ZPhiMCPWL")
chiExportFieldFunctionToVTKG(fflist[1],"ZPhiLBS")
