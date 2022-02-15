chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end

----############################################### Setup Transport mesh
tmesh = chiMeshHandlerCreate()

nodes={}
N=10
L=5.0
ds=L/N
xmin=0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end
mesh,region0 = chiMeshCreateUnpartitioned1DOrthoMesh(nodes)

chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)


tvol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,L-ds,L)


----############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,1.0,0.5)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

----############################################### Setup Transport Physics
chiMeshHandlerSetCurrent(tmesh)
phys0 = chiLBSCreateSolver()
chiSolverAddRegion(phys0,region0)

--========== Groups
grp = {}
for g=1,num_groups do
    grp[g] = chiLBSCreateGroup(phys0)
end

--========== ProdQuad
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,4)

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys0)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys0,gs0,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys0,gs0,pquad)
chiLBSGroupsetSetAngleAggDiv(phys0,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys0,cur_gs,1)
chiLBSGroupsetSetIterativeMethod(phys0,cur_gs,NPT_GMRES)
chiLBSGroupsetSetResidualTolerance(phys0,cur_gs,1.0e-6)
chiLBSGroupsetSetMaxIterations(phys0,cur_gs,300)
chiLBSGroupsetSetGMRESRestartIntvl(phys0,cur_gs,100)
--chiLBSGroupsetSetWGDSA(phys0,cur_gs,30,1.0e-4,false," ")
--chiLBSGroupsetSetTGDSA(phys0,cur_gs,30,1.0e-4,false," ")

--========== Boundary conditions
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
chiLBSSetProperty(phys0,BOUNDARY_CONDITION,
        ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys0,DISCRETIZATION_METHOD,PWLD)
chiLBSSetProperty(phys0,SCATTERING_ORDER,0)

chiLBSInitialize(phys0)
chiLBSExecute(phys0)
--
lbs_pwl_ff = chiGetFieldFunctionHandleByName("Flux_g0_m0")

----############################################### Setup Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys1 = chiMonteCarlonCreateSolver(chi_location_id, "FMCParTra")
chiSolverAddRegion(phys1,region0)

bsrc[1] = 1.0/4/math.pi
chiMonteCarlonCreateSource(phys1,"BOUNDARY_SOURCE",5,bsrc);

if (fac==nil) then
    fac=1
end
fv_offset = 0
fv_offset = num_groups
chiMonteCarlonSetProperty2(phys1,"NUM_PARTICLES"              ,fac*1e6)
chiMonteCarlonSetProperty2(phys1,"TALLY_MERGE_INTVL"          ,100e3)
chiMonteCarlonSetProperty2(phys1,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys1,"MONOENERGETIC"              ,true)
chiMonteCarlonSetProperty2(phys1,"FORCE_ISOTROPIC"            ,false)
chiMonteCarlonSetProperty2(phys1,"TALLY_MULTIPLICATION_FACTOR",1.0)
chiMonteCarlonSetProperty2(phys1,"MAKE_PWLD_SOLUTION"         ,true)
chiMonteCarlonSetProperty2(phys1,"PRINT_TFCS"                 ,true)

chiMonteCarlonAddCustomVolumeTally(phys1,tvol0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fmc_pwl_ff = chiGetFieldFunctionHandleByName("FMCParTra-PWLFlux_g0_m0")

----############################################### Setup Residual Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys2 = chiMonteCarlonCreateSolver(chi_location_id, "RMCParTra")
chiSolverAddRegion(phys2,region0)

bndry_sources = {{["bndry_index"] = ZMIN, ["bndry_src"] = {1}}}
chiMonteCarlonCreateSource(phys2,"RESIDUAL_TYPE_A",lbs_pwl_ff, bndry_sources);

chiMonteCarlonSetProperty2(phys2,"NUM_PARTICLES"              ,fac*1e6)
chiMonteCarlonSetProperty2(phys2,"TALLY_MERGE_INTVL"          ,100e3)
chiMonteCarlonSetProperty2(phys2,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys2,"MONOENERGETIC"              ,true)
chiMonteCarlonSetProperty2(phys2,"FORCE_ISOTROPIC"            ,true)
chiMonteCarlonSetProperty2(phys2,"TALLY_MULTIPLICATION_FACTOR",1.0/1.0)
chiMonteCarlonSetProperty2(phys2,"MAKE_PWLD_SOLUTION"         ,true)
chiMonteCarlonSetProperty2(phys2,"PRINT_TFCS"                 ,true)

chiMonteCarlonSetProperty2(phys2,"RESIDUAL_SRC_FF_OPTION"     ,2)
chiMonteCarlonSetProperty2(phys2,"RESIDUAL_SRC_NY"            ,100000)

chiMonteCarlonAddCustomVolumeTally(phys2,tvol0)

chiSolverInitialize(phys2)
chiSolverExecute(phys2)

rmc_pwl_ff = chiGetFieldFunctionHandleByName("RMCParTra-PWLFlux_g0_m0")

----############################################### Getting Sn and MC solution
cline0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline0,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline0,LINE_SECONDPOINT,0.0,0.0, L+xmin)
chiFFInterpolationSetProperty(cline0,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline0,ADD_FIELDFUNCTION,lbs_pwl_ff)
chiFFInterpolationSetProperty(cline0,ADD_FIELDFUNCTION,fmc_pwl_ff)

chiFFInterpolationInitialize(cline0)
chiFFInterpolationExecute(cline0)
chiFFInterpolationExportPython(cline0)

--Get values
arrays = chiFFInterpolationGetValue(cline0)
num_blocks = rawlen(arrays)

if (num_blocks>=2) then
    num_vals = rawlen(arrays[1])
    print("Number of values in first block: "..tostring(num_vals))

    MC_sol = arrays[2]
    SN_sol = arrays[1]

    true_error = {}
    for i=1,num_vals do
        true_error[i] = MC_sol[i] - SN_sol[i]
    end
end

----############################################### Getting RMC solution
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, L+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,rmc_pwl_ff)

chiFFInterpolationSetProperty(cline,LINE_CUSTOM_ARRAY,true_error)

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

----############################################### Show plots
if (chi_location_id == 0) then
    local handle = io.popen("python3 ZLFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")
end
