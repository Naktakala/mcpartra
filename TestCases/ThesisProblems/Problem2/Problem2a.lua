chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end

----############################################### Setup Transport mesh
tmesh = chiMeshHandlerCreate()

nodes={}
N=60
L=10.0
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

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,1.5,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

vol2 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,10*58/60,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol2,2)


tvol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,L-2*ds,L)


----############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material");
materials[3] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[3],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[3],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,1.0,0.0)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.2,0.2)
chiPhysicsMaterialSetProperty(materials[3],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.01,0.0)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 10.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[3],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE,1)

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
--bsrc={}
--for g=1,num_groups do
--    bsrc[g] = 0.0
--end
--bsrc[1] = 1.0/2
--chiLBSSetProperty(phys0,BOUNDARY_CONDITION,
--        ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

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

chiMonteCarlonCreateSource(phys1,"MATERIAL_SRC");

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
chiMonteCarlonSetProperty2(phys1,"PRINT_TFCS"                 ,false)
if (no_tr == nil) then no_tr = false; end
chiMonteCarlonSetProperty2(phys1,"NO_TRANSPORT"               ,no_tr)

chiMonteCarlonAddCustomVolumeTally(phys1,tvol0)
chiMonteCarlonAddCustomVolumeTally(phys1,vol0)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fmc_pwl_ff = chiGetFieldFunctionHandleByName("FMCParTra-PWLFlux_g0_m0")

----############################################### Setup Residual Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys2 = chiMonteCarlonCreateSolver(chi_location_id, "RMCParTra")
chiSolverAddRegion(phys2,region0)

bndry_sources = {}
chiMonteCarlonCreateSource(phys2,"RESIDUAL_TYPE_A",lbs_pwl_ff, bndry_sources);

chiMonteCarlonSetProperty2(phys2,"NUM_PARTICLES"              ,fac*1e6)
chiMonteCarlonSetProperty2(phys2,"TALLY_MERGE_INTVL"          ,100e3)
chiMonteCarlonSetProperty2(phys2,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys2,"MONOENERGETIC"              ,true)
chiMonteCarlonSetProperty2(phys2,"FORCE_ISOTROPIC"            ,true)
chiMonteCarlonSetProperty2(phys2,"TALLY_MULTIPLICATION_FACTOR",1.0/1.0)
chiMonteCarlonSetProperty2(phys2,"MAKE_PWLD_SOLUTION"         ,true)
chiMonteCarlonSetProperty2(phys2,"PRINT_TFCS"                 ,false)
if (no_tr == nil) then no_tr = false; end
chiMonteCarlonSetProperty2(phys2,"NO_TRANSPORT"               ,no_tr)

-- 0 PWLD
-- 1 Q0
-- 2 PWLC
-- 3 ZERO
if (r_opt == nil) then r_opt = 0; end
chiMonteCarlonSetProperty2(phys2,"RESIDUAL_SRC_FF_OPTION"     ,r_opt)
chiMonteCarlonSetProperty2(phys2,"RESIDUAL_SRC_NY"            ,100000)

chiMonteCarlonAddCustomVolumeTally(phys2,tvol0)
chiMonteCarlonAddCustomVolumeTally(phys2,vol0)

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

-- 0 PWLD
-- 1 Q0
-- 2 PWLC
-- 3 ZERO

--SMC     [0]  m=0 g=0 : 2.7154e-03 7.7338e-06 2.8481e-03  16420  1927
--phi_0   [0]  m=0 g=0 : 2.7052e-03 2.5564e-06 9.4499e-04   3943  1344

--Q0_S2   [0]  m=0 g=0 : 1.1795e-03 9.3239e-06 7.9047e-03   3047  1296
--PWLD_S2 [0]  m=0 g=0 : 1.0997e-03 5.5369e-06 5.0349e-03   2627  1280
--PWLC_S2 [0]  m=0 g=0 : 1.1117e-03 5.5342e-06 4.9780e-03   2605  1289

--Q0_S8   [0]  m=0 g=0 :  4.2639e-05 9.5021e-06  2.2285e-01  2871  1297
--PWLD_S8 [0]  m=0 g=0 : -1.3398e-05 5.8525e-06 -4.3683e-01  2595  1296
--PWLC_S8 [0]  m=0 g=0 :  2.7275e-05 5.7766e-06  2.1179e-01  2579  1194


