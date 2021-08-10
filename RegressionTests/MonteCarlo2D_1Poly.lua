-- 2D Transport test with Vacuum and Incident-isotropic BC.
-- SDM: PWLD
-- Test: Max-value=3.33722 and 2.49672e-01
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
chiMeshCreateUnpartitioned2DOrthoMesh(mesh,mesh)
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

chiMonteCarlonCreateSource(phys1,MCSrcTypes.MATERIAL_SRC);

chiMonteCarlonSetProperty2(phys1,"NUM_PARTICLES"              ,10e6)
chiMonteCarlonSetProperty2(phys1,"TALLY_MERGE_INTVL"          ,1e5)
chiMonteCarlonSetProperty2(phys1,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys1,"MONOENERGETIC"              ,false)
chiMonteCarlonSetProperty2(phys1,"FORCE_ISOTROPIC"            ,false)
chiMonteCarlonSetProperty2(phys1,"TALLY_MULTIPLICATION_FACTOR",1.0)
chiMonteCarlonSetProperty2(phys1,"MAKE_PWLD_SOLUTION"         ,true)
chiMonteCarlonSetProperty2(phys1,"UNCOLLIDED_ONLY"            ,false)

chiMonteCarlonInitialize(phys1)
chiMonteCarlonExecute(phys1)

mc_pwl_ff = chiGetFieldFunctionHandleByName("MCParTra-PWLFlux_g0_m0")
mc_pwl_ff_g1 = chiGetFieldFunctionHandleByName("MCParTra-PWLFlux_g1_m0")

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

--############################################### Get field functions
lbs_pwl_ff = chiGetFieldFunctionHandleByName("Flux_g0_m0")
lbs_pwl_ff_g1 = chiGetFieldFunctionHandleByName("Flux_g1_m0")

--############################################### Line plot
--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,-0.01,-L/2, 0.0)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,-0.01,L/2, 0.0)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,mc_pwl_ff)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,lbs_pwl_ff)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,mc_pwl_ff_g1)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,lbs_pwl_ff_g1)

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

--############################################### Volume integrations
ffi1 = chiFFInterpolationCreate(VOLUME)
curffi = ffi1
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol0)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,mc_pwl_ff)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value1=%.5f", maxval))

ffi2 = chiFFInterpolationCreate(VOLUME)
curffi = ffi2
chiFFInterpolationSetProperty(curffi,OPERATION,OP_MAX)
chiFFInterpolationSetProperty(curffi,LOGICAL_VOLUME,vol1)
chiFFInterpolationSetProperty(curffi,ADD_FIELDFUNCTION,mc_pwl_ff_g1)

chiFFInterpolationInitialize(curffi)
chiFFInterpolationExecute(curffi)
maxval = chiFFInterpolationGetValue(curffi)

chiLog(LOG_0,string.format("Max-value2=%.5e", maxval))

--############################################### Exports
chiExportFieldFunctionToVTKG(mc_pwl_ff,"ZPhiMCPWL")
chiExportFieldFunctionToVTKG(lbs_pwl_ff,"ZPhiLBS")

--############################################### Plots
if (chi_location_id == 0 and master_export==nil) then
    local handle = io.popen("python ZLFFI00.py")
end
