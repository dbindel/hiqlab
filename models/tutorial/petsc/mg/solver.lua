-- Loading the Lua input file
require 'meshes1.lua'

PetscViewerSetFormat(PETSC_VIEWER_STDOUT_WORLD,'matlab')
mypid = MPI_Comm_rank(PETSC_COMM_WORLD);
wsize = MPI_Comm_size(PETSC_COMM_WORLD);

-- Exit if not correct number of process
if (wsize~=nparts) then
    error('World size is not equal to partition size')
end

-- Set up to read partition
mnpr  = Mesh_Partition_Reader:new()

-- Construct meshes(numbered coarse to fine)
meshes={}
for i = 1,nlevels do
    meshes[i] = Mesh_Partition:new(mypid,ndm)
    dense_x = dense_xa[i]
    dense_y = dense_ya[i]
    meshes[i]:set_partition(mnpr,fnamea[i],mesh_func)
end

-- Assemble prolongators
mat_type = 2
is_reduced=1
P = {}
for i = 1,nlevels-1 do
    P[i] = Mesh_Partition_assemble_P_petsc2(meshes[i],meshes[i+1],fpnamea[i],mat_type,is_reduced)
end

-- Assemble stiffness
Kz = Mesh_Partition_assemble_dR_petsc(meshes[nlevels],1,0,0,mat_type,is_reduced)

-- Assemble the RHS
Fz = Mesh_Partition_assemble_R_petsc(meshes[nlevels],is_reduced)
VecScale(Fz,-1)

-- Assemble LHS
Uz = MatGetVecX(Kz)
VecSetFromOptions(Uz)

-- Create PC
pc = PCCreate()
PCSetOperators(pc,Kz,Kz,"DIFFERENT_NONZERO_PATTERN");
PCSetType(pc,"mg")
PCMGSetLevels(pc,nlevels)
PCMGSetType(pc,"PC_MG_MULTIPLICATIVE")
PCMGSetCycleType(pc,"PC_MG_CYCLE_V");
PCMGSetInterpolates(pc,P)
PCRegisterKaczmarz()
--PCMGSetStationarySmoothers(pc,"kaczmarz",1.0,5,1, 1, 1)
PCMGSetStationarySmoothers(pc,"sor",1.0,1,1,5,5)
--PCMGSetKrylovSmoothers(pc,"gmres","none",5,5,1.0,10,1,1)
--PCMGSetKrylovSmoothers(pc,"gmres","kaczmarz",5,5,1.0,15,1,1)
--PCMGSetDirectCoarseSolve(pc,"superlu_dist")
PCMGSetDirectCoarseSolve(pc,"aijmumps")
PCView(pc,PETSC_VIEWER_STDOUT_WORLD)
--Acoarse = PCMGGetOperator(pc,0)
--petsc_matrix_dump(Acoarse,'Acoarse.m')

-- Create KSP
ksp = KSPCreate()
KSPSetOperators(ksp, Kz, Kz, "DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-9, PETSC_DEFAULT, PETSC_DEFAULT, 20)
PetscOptionsSetValue("-ksp_type","gmres")
KSPSetPC(ksp,pc)

-- Solve
KSPSetFromOptions(ksp)
KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD)
KSPSolve(ksp,Fz,Uz)

-- print results
pv = PetscViewerASCIIOpen('petscmat1.m')
PetscViewerSetFormat(pv,'matlab')
MatView(Kz,pv)
VecView(Uz,pv)
VecView(Fz,pv)
PetscViewerDestroy(pv)

-- Results
rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)
VecScale(Fz,-1)
norm2 = qVecNorm(Fz,2)
MatMultAdd(Kz,Uz,Fz,Fz)
norm1 = qVecNorm(Fz,2)
if (mypid==0) then
    print('Residual:',rnorm, ' NumIterations:', its)
    print('Residual(ex):', norm1/norm2, ' NumIterations:', its)
end

-- clean up
KSPDestroy(ksp)
PCDestroy(pc)
VecDestroy(Fz)
VecDestroy(Uz)
for i = 1,nlevels-1 do
    MatDestroy(P[i])
end
MatDestroy(Kz)
for i = 1,nlevels do
    meshes[i]:delete()
end
mnpr:delete()
