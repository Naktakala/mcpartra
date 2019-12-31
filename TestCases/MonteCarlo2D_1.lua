--G=168
chiMPIBarrier()
if (chi_location_id == 0) then
    print("############################################### LuaTest")
end
--dofile(CHI_LIBRARY)



--############################################### Setup mesh
chiMeshHandlerCreate()

newSurfMesh = chiSurfaceMeshCreate();
chiSurfaceMeshImportFromOBJFile(newSurfMesh,
        "CHI_RESOURCES/TestObjects/SquareMesh2x2Quads.obj",true)

--############################################### Extract edges from surface mesh
loops,loop_count = chiSurfaceMeshGetEdgeLoopsPoly(newSurfMesh)

line_mesh = {};
line_mesh_count = 0;

for k=1,loop_count do
    split_loops,split_count = chiEdgeLoopSplitByAngle(loops,k-1);
    for m=1,split_count do
        line_mesh_count = line_mesh_count + 1;
        line_mesh[line_mesh_count] =
        chiLineMeshCreateFromLoop(split_loops,m-1);
    end

end

--############################################### Setup Regions
region1 = chiRegionCreate()
chiRegionAddSurfaceBoundary(region1,newSurfMesh);
for k=1,line_mesh_count do
    chiRegionAddLineBoundary(region1,line_mesh[k]);
end

--############################################### Create meshers
chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED);
chiVolumeMesherCreate(VOLUMEMESHER_PREDEFINED2D);

chiVolumeMesherSetProperty(FORCE_POLYGONS,true);
--chiVolumeMesherSetProperty(MESH_GLOBAL,true)

chiSurfaceMesherSetProperty(PARTITION_X,2)
chiSurfaceMesherSetProperty(PARTITION_Y,2)
chiSurfaceMesherSetProperty(CUT_X,0.0)
chiSurfaceMesherSetProperty(CUT_Y,0.0)


--############################################### Execute meshing
chiSurfaceMesherExecute();
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


num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_3_170.data")
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        PDT_XSFILE,"CHI_TEST/xs_3_170.data")

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS0,num_groups,0.1)


src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)



--############################################### Setup Physics
phys1 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys1,region1)

chiMonteCarlonCreateSource(phys1,MCSrcTypes.BNDRY_SRC,1);

chiMonteCarlonSetProperty(phys1,MCProperties.NUM_PARTICLES,10e6)
chiMonteCarlonSetProperty(phys1,MCProperties.TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MERGE_INTVL,1e5)
chiMonteCarlonSetProperty(phys1,MCProperties.SCATTERING_ORDER,10)
chiMonteCarlonSetProperty(phys1,MCProperties.MONOENERGETIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.FORCE_ISOTROPIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0)
chiMonteCarlonSetProperty(phys1,MCProperties.MAKE_PWLD_SOLUTION,true)
chiMonteCarlonSetProperty(phys1,MCProperties.UNCOLLIDED_ONLY,false)

chiMonteCarlonInitialize(phys1)
chiMonteCarlonExecute(phys1)


--Testing consolidated interpolation
cline = chiFFInterpolationCreate(LINE)
chiFFInterpolationSetProperty(cline,LINE_FIRSTPOINT,0.0,-1.0, 0.0)
chiFFInterpolationSetProperty(cline,LINE_SECONDPOINT,0.0,1.0, 0.0)
chiFFInterpolationSetProperty(cline,LINE_NUMBEROFPOINTS, 500)

for k=1,2 do
    chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,k-1)
end
chiFFInterpolationSetProperty(cline,ADD_FIELDFUNCTION,168)


chiFFInterpolationInitialize(cline)
chiFFInterpolationExecute(cline)
chiFFInterpolationExportPython(cline)


--


if (chi_location_id == 0) then
    local handle = io.popen("python ZLFFI00.py")
end

chiExportFieldFunctionToVTKG(168,"ZPhiMC")