chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end

----############################################### Setup Transport mesh
tmesh = chiMeshHandlerCreate()

nodes={}
N=60
L=5.0
ds=L/N
xmin=0.0
for i=0,N do
    nodes[i+1] = xmin + i*ds
end
mesh,region0 = chiMeshCreateUnpartitioned2DOrthoMesh(nodes,nodes)

--chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_XYZ_STYLE)
chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)

--Execute meshing
if (chi_number_of_processes == 4) then
    chiSurfaceMesherSetProperty(PARTITION_X,2)
    chiSurfaceMesherSetProperty(PARTITION_Y,2)
    chiSurfaceMesherSetProperty(CUT_X,L/2)
    chiSurfaceMesherSetProperty(CUT_Y,L/2)
end
chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,L/3,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)





----############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,1.0,0.0)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,1.0,0.0)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 3.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)






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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,1,1)

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

if (fac == nil) then fac=0.01 end
fv_offset = 0
fv_offset = num_groups

chiMonteCarlonSetProperty2(phys1,"NUM_PARTICLES"              ,fac*100e6)
chiMonteCarlonSetProperty2(phys1,"TALLY_MERGE_INTVL"          ,100e3)
chiMonteCarlonSetProperty2(phys1,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys1,"MONOENERGETIC"              ,true)
chiMonteCarlonSetProperty2(phys1,"FORCE_ISOTROPIC"            ,false)
chiMonteCarlonSetProperty2(phys1,"TALLY_MULTIPLICATION_FACTOR",1.0)
chiMonteCarlonSetProperty2(phys1,"MAKE_PWLD_SOLUTION"         ,true)

xmin=2.33333;ymin=4.16666;dx=0.33333;dy=0.16666
tvol0 = chiLogicalVolumeCreate(RPP,xmin,xmin+dx,ymin,ymin+dy,-1000,1000)
xmin=2.33333;ymin=2.33333;dx=0.33333;dy=0.16666
tvol1 = chiLogicalVolumeCreate(RPP,xmin,xmin+dx,ymin,ymin+dy,-1000,1000)

--tvol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)

chiMonteCarlonAddCustomVolumeTally(phys1,tvol0)
chiMonteCarlonAddCustomVolumeTally(phys1,tvol1)

chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fmc_pwl_ff = chiGetFieldFunctionHandleByName("FMCParTra-PWLFlux_g0_m0")

----############################################### Setup ref Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys2 = chiMonteCarlonCreateSolver(chi_location_id, "RMCParTra")
chiSolverAddRegion(phys2,region0)

chiMonteCarlonCreateSource(phys2,"RESIDUAL_TYPE_A",lbs_pwl_ff);

chiMonteCarlonSetProperty2(phys2,"NUM_PARTICLES"              ,fac*100e6)
chiMonteCarlonSetProperty2(phys2,"TALLY_MERGE_INTVL"          ,100e3)
chiMonteCarlonSetProperty2(phys2,"SCATTERING_ORDER"           ,0)
chiMonteCarlonSetProperty2(phys2,"MONOENERGETIC"              ,true)
chiMonteCarlonSetProperty2(phys2,"FORCE_ISOTROPIC"            ,false)
chiMonteCarlonSetProperty2(phys2,"TALLY_MULTIPLICATION_FACTOR",1.0)
chiMonteCarlonSetProperty2(phys2,"MAKE_PWLD_SOLUTION"         ,true)

chiMonteCarlonAddCustomVolumeTally(phys2,tvol0)
chiMonteCarlonAddCustomVolumeTally(phys2,tvol1)

chiSolverInitialize(phys2)
chiSolverExecute(phys2)

rmc_pwl_ff = chiGetFieldFunctionHandleByName("RMCParTra-PWLFlux_g0_m0")

----############################################### Getting Sn and MC solution
cline0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline0,LINE_FIRSTPOINT ,L/2+0.001,0.0,0.0)
chiFFInterpolationSetProperty(cline0,LINE_SECONDPOINT,L/2+0.001,0.0+L,0.0)
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
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT ,L/2,0.0,0.0)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,L/2,0.0+L,0.0)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,rmc_pwl_ff)

chiFFInterpolationSetProperty(cline,LINE_CUSTOM_ARRAY,true_error)

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

----############################################### Slice
slice2 = chiFFInterpolationCreate(SLICE)
chiFFInterpolationSetProperty(slice2,SLICE_POINT,0.0,0.0,0.025)
chiFFInterpolationSetProperty(slice2,ADD_FIELDFUNCTION,rmc_pwl_ff)

chiFFInterpolationInitialize(slice2)
chiFFInterpolationExecute(slice2)
chiFFInterpolationExportPython(slice2)

----############################################### Show plots
if ((chi_location_id == 0) and (with_plot~=nil)) then
    local handle = io.popen("python3 ZLFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")
    local handle = io.popen("python3 ZPFFI20.py")
end

chiExportFieldFunctionToVTKG(rmc_pwl_ff,"ZRMC")
chiExportFieldFunctionToVTKG(fmc_pwl_ff,"ZFMC")
chiExportFieldFunctionToVTKG(lbs_pwl_ff,"ZSn")
