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
mesh,region0 = chiMeshCreate1DSlabMesh(nodes)

chiVolumeMesherSetProperty(PARTITION_Z,chi_number_of_processes)

--Execute meshing
chiVolumeMesherExecute();

----############################################### Set Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,L/3,1000)
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
bsrc={}
for g=1,num_groups do
    bsrc[g] = 0.0
end
bsrc[1] = 1.0/2
--chiLBSSetProperty(phys0,BOUNDARY_CONDITION,
--        ZMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);

--========== Solvers
chiLBSSetProperty(phys0,PARTITION_METHOD,FROM_SURFACE)
chiLBSSetProperty(phys0,DISCRETIZATION_METHOD,PWLD3D)
chiLBSSetProperty(phys0,SCATTERING_ORDER,0)

chiLBSInitialize(phys0)
chiLBSExecute(phys0)
--
fflist0,count = chiLBSGetScalarFieldFunctionList(phys0)


----############################################### Setup Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys1 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys1,region0)

chiMonteCarlonCreateSource(phys1,MCSrcTypes.MATERIAL_SRC);

fac=1
fv_offset = 0
fv_offset = num_groups
chiMonteCarlonSetProperty(phys1,MCProperties.NUM_PARTICLES,fac*0.5e6)
chiMonteCarlonSetProperty(phys1,MCProperties.TFC_UPDATE_INTVL,100e3)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MERGE_INTVL,100e3)
chiMonteCarlonSetProperty(phys1,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys1,MCProperties.MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys1,MCProperties.FORCE_ISOTROPIC,false)
--chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MULTIPLICATION_FACTOR,5.0*3/3)
chiMonteCarlonSetProperty(phys1,MCProperties.MAKE_PWLD_SOLUTION,true)

chiMonteCarlonInitialize(phys1)
chiMonteCarlonExecute(phys1)

fflist1,count = chiGetFieldFunctionList(phys1)

----############################################### Setup ref Monte Carlo Physics
chiMeshHandlerSetCurrent(tmesh)
phys2 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys2,region0)

chiMonteCarlonCreateSource(phys2,MCSrcTypes.RESIDUAL_TYPE_B,-1,fflist0[1],0.0);


chiMonteCarlonSetProperty(phys2,MCProperties.NUM_UNCOLLIDED_PARTICLES,fac*1e6)
chiMonteCarlonSetProperty(phys2,MCProperties.NUM_PARTICLES,fac*0.5*1e6)
chiMonteCarlonSetProperty(phys2,MCProperties.TFC_UPDATE_INTVL,100e3)
chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MERGE_INTVL,100e3)
chiMonteCarlonSetProperty(phys2,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys2,MCProperties.MONOENERGETIC,true)
chiMonteCarlonSetProperty(phys2,MCProperties.FORCE_ISOTROPIC,true)
--chiMonteCarlonSetProperty(phys2,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0/1.0)
chiMonteCarlonSetProperty(phys2,MCProperties.MAKE_PWLD_SOLUTION,true)


chiMonteCarlonInitialize(phys2)
chiMonteCarlonExecute(phys2)

fflist2,count = chiGetFieldFunctionList(phys2) --Fine mesh MC

----############################################### Getting Sn and MC solution
cline0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline0,LINE_FIRSTPOINT,0.0,0.0,0.0+xmin)
chiFFInterpolationSetProperty(cline0,LINE_SECONDPOINT,0.0,0.0, 5.0+xmin)
chiFFInterpolationSetProperty(cline0,LINE_NUMBEROFPOINTS, 500)

chiFFInterpolationSetProperty(cline0,ADD_FIELDFUNCTION,fflist0[1])
chiFFInterpolationSetProperty(cline0,ADD_FIELDFUNCTION,fflist1[1]+num_groups)

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
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,0.0, 5.0+xmin)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

--chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1]+num_groups)
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,fflist2[1])

chiFFInterpolationSetProperty(cline,LINE_CUSTOM_ARRAY,true_error)

chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)

----############################################### Show plots
if (chi_location_id == 0) then
    local handle = io.popen("python3 ZLFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")
end
