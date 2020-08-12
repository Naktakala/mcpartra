if (dim == nil) then dim = 2; end

--############################################### Mesh
chiMeshHandlerCreate()

if (L == nil) then L = 5.0; end
if (N == nil) then N = 40; end

ds = L/N

nodes={}
for i=0,N do
    nodes[i+1] = i*ds
end

if (dim == 2) then
    chiMeshCreate2DOrthoMesh(nodes,nodes)
    chiVolumeMesherExecute();

    --########################################## Material IDs
    vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
    chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

    vol1 = chiLogicalVolumeCreate(RPP,0.0,0.5*10/N,
                                      0.0,0.5*10/N,
                                      -1000,1000)
    chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)
end

if (dim == 3) then
    chiMeshCreate3DOrthoMesh(nodes,nodes,nodes)
    chiVolumeMesherExecute();

    --########################################## Material IDs
    vol0 = chiLogicalVolumeCreate(RPP,-1000,1000,-1000,1000,-1000,1000)
    chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol0,0)

    vol1 = chiLogicalVolumeCreate(RPP,0.0,0.5*10/N,
                                      0.0,0.5*10/N,
                                      0.0,0.5*10/N)
    chiVolumeMesherSetProperty(MATID_FROMLOGICAL,vol1,1)
end



--############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");
materials[2] = chiPhysicsAddMaterial("Test Material2");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)
chiPhysicsMaterialAddProperty(materials[2],TRANSPORT_XSECTIONS)

chiPhysicsMaterialAddProperty(materials[1],ISOTROPIC_MG_SOURCE)
chiPhysicsMaterialAddProperty(materials[2],ISOTROPIC_MG_SOURCE)


num_groups = 2
--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"CHI_TEST/xs_3_170.data")
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
--        PDT_XSFILE,"CHI_TEST/xs_3_170.data")

chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,
                                           SIMPLEXS0,num_groups,1.0)
chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,
                                           SIMPLEXS0,num_groups,1.0)

--chiPhysicsMaterialSetProperty(materials[1],TRANSPORT_XSECTIONS,SIMPLEXS1,num_groups,1.0,0.8)
--chiPhysicsMaterialSetProperty(materials[2],TRANSPORT_XSECTIONS,SIMPLEXS1,num_groups,1.0,0.8)

src={}
for g=1,num_groups do
    src[g] = 0.0
end
--src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[1],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)
src[1] = 1.0
chiPhysicsMaterialSetProperty(materials[2],ISOTROPIC_MG_SOURCE,FROM_ARRAY,src)


--############################################### Setup Physics
phys1 = chiMonteCarlonCreateSolver()
chiSolverAddRegion(phys1,region1)

--chiMonteCarlonCreateSource(phys1,MCSrcTypes.BNDRY_SRC,1);
chiMonteCarlonCreateSource(phys1,MCSrcTypes.MATERIAL_SRC);

chiMonteCarlonSetProperty(phys1,MCProperties.NUM_PARTICLES,20.0e6)
chiMonteCarlonSetProperty(phys1,MCProperties.TFC_UPDATE_INTVL,10e3)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MERGE_INTVL,1e5)
chiMonteCarlonSetProperty(phys1,MCProperties.SCATTERING_ORDER,0)
chiMonteCarlonSetProperty(phys1,MCProperties.MONOENERGETIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.FORCE_ISOTROPIC,false)
chiMonteCarlonSetProperty(phys1,MCProperties.TALLY_MULTIPLICATION_FACTOR,(0.5*10/N)^dim)
chiMonteCarlonSetProperty(phys1,MCProperties.MAKE_PWLD_SOLUTION,true)
chiMonteCarlonSetProperty(phys1,MCProperties.UNCOLLIDED_ONLY,true)

chiMonteCarlonInitialize(phys1)
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
pquad = chiCreateProductQuadrature(GAUSS_LEGENDRE_CHEBYSHEV,13, 64)
if (dim == 3) then
    --pquad = chiCreateSLDFESQAngularQuadrature(3)
end

--========== Groupset def
gs0 = chiLBSCreateGroupset(phys1)
cur_gs = gs0
chiLBSGroupsetAddGroups(phys1,cur_gs,0,num_groups-1)
chiLBSGroupsetSetQuadrature(phys1,cur_gs,pquad)
chiLBSGroupsetSetAngleAggDiv(phys1,cur_gs,1)
chiLBSGroupsetSetGroupSubsets(phys1,cur_gs,1)
chiLBSGroupsetSetAngleAggregationType(phys1,cur_gs,LBSGroupset.ANGLE_AGG_POLAR)
chiLBSGroupsetSetIterativeMethod(phys1,cur_gs,NPT_GMRES)
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
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        YMIN,LBSBoundaryTypes.INCIDENT_ISOTROPIC,bsrc);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        XMIN,LBSBoundaryTypes.REFLECTING);
--chiLBSSetProperty(phys1,BOUNDARY_CONDITION,
--        XMAX,LBSBoundaryTypes.REFLECTING);

--========== Solvers
chiLBSSetProperty(phys1,DISCRETIZATION_METHOD,PWLD2D)
chiLBSSetProperty(phys1,SCATTERING_ORDER,0)

chiLBSInitialize(phys1)
if (skipLBS == nil) then chiLBSExecute(phys1) end

fflist,count = chiLBSGetScalarFieldFunctionList(phys1)


chiExportFieldFunctionToVTKG(num_groups,"ZPhiMC")
if (skipLBS == nil) then
    chiExportFieldFunctionToVTKG(fflist[1],"ZPhiLBS")
end
