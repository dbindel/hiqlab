-- Loading the Lua input file
require 'mesh1.lua'

-- Specify parameters
order  =  1
dense_x=  2.0e0
dense_y=  1.0e-0

-- Run the Lua input file for single Mesh
mesh = Mesh:new(2);
mesh = mesh_func(mesh)
mesh:initialize()

-- assemble stiffness in PETSc
mat_type = 2
is_reduced=1
Kz = mesh:assemble_dR_petsc(1,0,0,mat_type,is_reduced)

-- assemble rhs in PETSc
Fz = mesh:assemble_R_petsc(is_reduced)

-- create direct solver
ksp = KSPCreate()
KSPSetType(ksp,'preonly')
pc  = KSPGetPC(ksp)
PCSetType(pc,'lu')

Kz2= MatConvert(Kz,'superlu_dist')
MatDestroy(Kz)
Kz = Kz2
PCSetOperators(pc,Kz,Kz,'DIFFERENT_NONZERO_PATTERN')

-- create solution vector and solve
Uz = MatGetVecX(Kz)
KSPSolve(ksp,Fz,Uz)
rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)
print('Residual:',rnorm, ' NumIterations:', its)

-- Print results
--pv = PETSC_VIEWER_STDOUT_WORLD
pv = PetscViewerASCIIOpen('petscmat.m')
PetscViewerSetFormat(pv,'matlab')
MatView(Kz,pv)
VecView(Fz,pv)
VecView(Uz,pv)
PetscViewerDestroy(pv)

-- clean up
KSPDestroy(ksp)
MatDestroy(Kz)
VecDestroy(Fz)
VecDestroy(Uz)

-- clean up
mesh:delete()
