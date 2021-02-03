print("Executing input file")
--############################################### Mesh
chiMeshHandlerCreate()

if (L == nil) then L = 5.0; end
if (N == nil) then N = 20; end
if (X == nil) then X = 0.0; end
if (s == nil) then s = 0; end

ds = L/N

nodes={}
for i=0,N do
    nodes[i+1] = i*ds
end

chiMeshCreateUnpartitioned3DOrthoMesh(nodes,nodes,nodes)
chiVolumeMesherExecute();

--########################################## Material IDs
vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

vol1 = chiLogicalVolumeCreate(RPP,0.0,0.5*10/N,
        0.0,0.5*10/N,
        0.0,0.5*10/N)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)

vol2 = chiLogicalVolumeCreate(RPP,0.0+X,1.0+X,
        0.5,1.5,
        0.0,1.0)
chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol2,2)

--############################################### Add materials
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
chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
        SIMPLEXS0,num_groups,1.0)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
        SIMPLEXS0,num_groups,1.0)
chiPhysicsMaterialSetProperty(materials[3],TRANSPORT_XSECTIONS,
        SIMPLEXS0,num_groups,50.0)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0/4.0/math.pi/2
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 0.0
chiPhysicsMaterialSetProperty(materials[3],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)

--########################################## Setup solver
phys0 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys0,region1)

chiMonteCarlonCreateSource(phys0,MCSrcTypes.MATERIAL_SRC);

chiMonteCarlonSetProperty(phys0,MCProperties.NUM_PARTICLES,5e6)
chiMonteCarlonSetProperty(phys0,MCProperties.TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MERGE_INTVL,1e5)
chiMonteCarlonSetProperty(phys0,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys0,MCProperties.MONOENERGETIC,false)
chiMonteCarlonSetProperty(phys0,MCProperties.FORCE_ISOTROPIC,false)
chiMonteCarlonSetProperty(phys0,MCProperties.TALLY_MULTIPLICATION_FACTOR,1.0)
chiMonteCarlonSetProperty(phys0,MCProperties.MAKE_PWLD_SOLUTION,true)
chiMonteCarlonSetProperty(phys0,MCProperties.UNCOLLIDED_ONLY,true)

chiMonteCarlonInitialize(phys0)
chiMonteCarlonExecute(phys0)

fflist,count = chiGetFieldFunctionList(phys0)
chiExportFieldFunctionToVTKG(fflist[1]+num_groups,"UPhiMC")