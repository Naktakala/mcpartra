chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end

----############################################### Setup Transport mesh
tmesh = chiMeshHandlerCreate()

xnodes={}
ynodes={}
Nx=60
Ny=120
Lx=5.0
Ly=10.0
dx=Lx/Nx
dy=Ly/Ny
xmin=0.0
ymin=-5.0
for i=0,Nx do
    xnodes[i+1] = xmin + i*dx
end
for i=0,Ny do
    ynodes[i+1] = ymin + i*dy
end
mesh,region0 = chiMeshCreateUnpartitioned2DOrthoMesh(xnodes,ynodes)

--chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_XYZ_STYLE)
chiVolumeMesherSetProperty(PARTITION_TYPE,KBA_STYLE_XYZ)

--Execute meshing
if (chi_number_of_processes == 4) then
    chiSurfaceMesherSetProperty(PARTITION_X,2)
    chiSurfaceMesherSetProperty(PARTITION_Y,2)
    chiSurfaceMesherSetProperty(CUT_X,Lx/2)
    chiSurfaceMesherSetProperty(CUT_Y,Ly/2)
end
chiVolumeMesherExecute();

----############################################### Set Material IDs
--All air
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

--Scattering medium y-extent
vol1 = chiLogicalVolumeCreate(RPP,-1000,1000,-0.4*Ly,0.4*Ly,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)



----############################################### Set Material IDs
----Air gap
--vol0b = chiLogicalVolumeCreate(RPP,-0.166666+2.5,0.166666+2.5,-1000,1000,-1000,1000)
--chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0b,0)
--
----Source
vol2 = chiLogicalVolumeCreate(RPP,-1000,1000,-0.166666*5,0.166666*5,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol2,2)

--End of gap scatterer
vol1b = chiLogicalVolumeCreate(RPP,-1+2.5,1+2.5,0.9*Ly/2,Ly/2,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1b,1)

vol1c = chiLogicalVolumeCreate(RPP,-1+2.5,1+2.5,-Ly/2,-0.9*Ly/2,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1c,1)

--chiRegionExportMeshToVTK(region0, "ZMeshAfter")
--os.exit()

----############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");
materials[3] = chiPhysicsAddMaterial("Test Material3");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[3],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[3],ISOTROPIC_MG_SOURCE)


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.01,0.01)
chiPhysicsMaterialSetProperty(materials[2],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.1*20,0.8)
chiPhysicsMaterialSetProperty(materials[3],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.1*20,0.8)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 3.0
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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,6,6)

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
chiMonteCarlonSetProperty2(phys1,"APPLY_SOURCE_IMPORTANCE_SAMPLING",A or B)
chiMonteCarlonSetProperty2(phys1,"APPLY_SOURCE_ANGULAR_BIASING",B)

tvol0 = chiLogicalVolumeCreate(RPP,2.3333,2.6666,4.16666,4.33333,-1000,1000)
tvol1 = chiLogicalVolumeCreate(RPP,0.5   ,0.8333,4.16666,4.33333,-1000,1000)

--tvol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
--tvol1 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)

fmc_QOI_left = chiMonteCarlonAddCustomVolumeTally(phys1,tvol0)
fmc_QOI_rite = chiMonteCarlonAddCustomVolumeTally(phys1,tvol1)

chiMonteCarlonReadImportanceMap(phys1, "/Users/janv4/Desktop/ChiTech/LBAdjointSolver/Residual3MC2D_2b.o")
chiSolverInitialize(phys1)
chiSolverExecute(phys1)

fmc_pwl_ff = chiGetFieldFunctionHandleByName("FMCParTra-PWLFlux_g0_m0")

fmc_QOI_left_val, fmc_QOI_left_sig = chiMonteCarlonGetCustomVolumeTallyValue(phys1, fmc_QOI_left, 0, 0)
fmc_QOI_rite_val, fmc_QOI_rite_sig = chiMonteCarlonGetCustomVolumeTallyValue(phys1, fmc_QOI_rite, 0, 0)

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
chiMonteCarlonSetProperty2(phys2,"APPLY_SOURCE_IMPORTANCE_SAMPLING",A or B)
chiMonteCarlonSetProperty2(phys2,"APPLY_SOURCE_ANGULAR_BIASING",B)

rmc_QOI_left = chiMonteCarlonAddCustomVolumeTally(phys2,tvol0)
rmc_QOI_rite = chiMonteCarlonAddCustomVolumeTally(phys2,tvol1)

chiMonteCarlonReadImportanceMap(phys2, "/Users/janv4/Desktop/ChiTech/LBAdjointSolver/Residual3MC2D_2b.o")
chiSolverInitialize(phys2)
chiMonteCarlonExportImportanceMap(phys2, "Z_")
chiSolverExecute(phys2)
-- chiSolverExecute(phys1)

rmc_pwl_ff = chiGetFieldFunctionHandleByName("RMCParTra-PWLFlux_g0_m0")

rmc_QOI_left_val, rmc_QOI_left_sig = chiMonteCarlonGetCustomVolumeTallyValue(phys2, rmc_QOI_left, 0, 0)
rmc_QOI_rite_val, rmc_QOI_rite_sig = chiMonteCarlonGetCustomVolumeTallyValue(phys2, rmc_QOI_rite, 0, 0)

----############################################### Getting Sn and MC solution
cline0 = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline0,LINE_FIRSTPOINT ,Lx/2+0.001,ymin,0.0)
chiFFInterpolationSetProperty(cline0,LINE_SECONDPOINT,Lx/2+0.001,ymin+Ly,0.0)
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
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT ,Lx/2,ymin,0.0)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,Lx/2,ymin+Ly,0.0)
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
print(with_plot)
if ((chi_location_id == 0) and (with_plot==nil)) then
    local handle = io.popen("python3 ZLFFI00.py")
    local handle = io.popen("python3 ZLFFI10.py")
    local handle = io.popen("python3 ZPFFI20.py")
end

chiExportFieldFunctionToVTKG(rmc_pwl_ff,"ZRMC")
chiExportFieldFunctionToVTKG(fmc_pwl_ff,"ZFMC")
chiExportFieldFunctionToVTKG(lbs_pwl_ff,"ZSn")

print(" ")
print("FMC-QOI-left "..string.format("%.4e",fmc_QOI_left_val).." "..string.format("%.4e",fmc_QOI_left_sig))
print("FMC-QOI-rite "..string.format("%.4e",fmc_QOI_rite_val).." "..string.format("%.4e",fmc_QOI_rite_sig))
print(" ")
print("RMC-QOI-left "..string.format("%.4e",rmc_QOI_left_val).." "..string.format("%.4e",rmc_QOI_left_sig))
print("RMC-QOI-rite "..string.format("%.4e",rmc_QOI_rite_val).." "..string.format("%.4e",rmc_QOI_rite_sig))
print(" ")

--
--FMC-QOI-left 1.1142e-02 1.0899e-04
--FMC-QOI-rite 9.0198e-04 3.1225e-05
--
--RMC-QOI-left 6.3262e-03 4.5375e-05
--RMC-QOI-rite 1.1938e-04 1.9710e-05

--Cell-importance, No angular-importance
--FMC-QOI-left 1.1423e-02 8.9558e-05 (1.48)
--FMC-QOI-rite 8.6983e-04 2.4800e-05 (1.59)
--
--RMC-QOI-left 7.0237e-03 3.6239e-05 (1.57)
--RMC-QOI-rite 3.9783e-04 1.2472e-05 (2.50)

--Cell-importance, And angular-importance
--FMC-QOI-left 1.1140e-02 6.2713e-05 (3.02)
--FMC-QOI-rite 8.9446e-04 1.7995e-05 (3.01)
--
--RMC-QOI-left 7.4333e-03 3.6122e-05 (1.58)
--RMC-QOI-rite 5.1501e-04 1.0178e-05 (3.75)
