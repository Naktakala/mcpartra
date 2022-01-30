chiMeshHandlerCreate()
umesh = chiCreateEmptyUnpartitionedMesh()

chiUnpartitionedMeshUploadVertex(umesh, 0-0.5, 0-0.5, 0-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 1-0.5, 0-0.5, 0-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 1-0.5, 1-0.5, 0-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 0-0.5, 1-0.5, 0-0.5)

chiUnpartitionedMeshUploadVertex(umesh, 0-0.5, 0-0.5, 1-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 1-0.5, 0-0.5, 1-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 1-0.5, 1-0.5, 1-0.5)
chiUnpartitionedMeshUploadVertex(umesh, 0-0.5, 1-0.5, 1-0.5)

cell = {}
cell.type        = "POLYHEDRON"
cell.sub_type    = "HEXAHEDRON"
cell.num_faces   = 6
cell.material_id = 0
cell.face0 = {1,2,6,5}
cell.face1 = {0,4,7,3}
cell.face2 = {2,3,7,6}
cell.face3 = {0,1,5,4}
cell.face4 = {4,5,6,7}
cell.face5 = {0,3,2,1}

chiUnpartitionedMeshUploadCell(umesh, cell, true)
chiUnpartitionedMeshFinalizeEmpty(umesh)

region1 = chiRegionCreate()

chiSurfaceMesherCreate(SURFACEMESHER_PREDEFINED)
chiVolumeMesherCreate(VOLUMEMESHER_UNPARTITIONED)

chiSurfaceMesherExecute()
chiVolumeMesherExecute()

----############################################### Add materials
materials = {}
materials[1] = chiPhysicsAddMaterial("Test Material");

chiPhysicsMaterialAddProperty(materials[1],TRANSPORT_XSECTIONS)

num_groups = 1
chiPhysicsMaterialSetProperty(materials[1],
                              TRANSPORT_XSECTIONS,
                              SIMPLEXS1,1,0.01,0.01)


--#########################################################
if (seed == nil) then seed=0; end
phys1 = chiMonteCarlonCreateSolver(seed, "FMCParTra")
chiSolverAddRegion(phys1,region1)

chiSolverInitialize(phys1)
--chiSolverExecute(phys1)

TestCode(phys1);