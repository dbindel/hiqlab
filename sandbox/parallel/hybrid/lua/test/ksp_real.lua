ts   = os.clock()
require '../../lua_petsc/petsc-util.lua'

-- Loading the Lua input file
f0    = 0
order = 1;
dense = 1e-6; -- 1e-6
--runfile = loadfile('../../lua_trilinos/test/block3d.lua')
--runfile = loadfile('pml1d.lua')
runfile = loadfile('plate3d.lua')
--runfile = loadfile('/home/tkoyama/programs/hiqlab/models/developing/disk_resonator/3d/disk3d.lua')


-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
---[[
mesh2 = mesh:get_clean_mesh()
mesh2:set_bc(bc_function)
mesh2:initialize()
mesh2:apply_bc()
--]]
numid = mesh2:get_numid()
print("Numid:",numid)

-- Assemble distributed stiffness
mattype    = 3;
is_reduced = 2;
w0         = 0
--w0         = 15991068*2*pi
Kz = mesh2:assemble_dR_petsc(1, 0, -w0*w0, mattype, is_reduced)

-- Assemble the RHS
Fz = mesh2:assemble_R_petsc(is_reduced)
VecScale(Fz,-1)

-- Assemble LHS
Uz = MatGetVecX(Kz)
VecSetFromOptions(Uz)

-- Create KSP
ksp = KSPCreate()
KSPSetOperators(ksp, Kz, Kz, "DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-7, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT)
PetscOptionsSetValue("-ksp_type","gmres")
PetscOptionsSetValue("-ksp_monitor","")
--PetscOptionsSetValue("-log_summary","")
pc  = KSPGetPC(ksp)

-- HYPRE PC
--[[
PCSetType(pc, "hypre")
PCHYPRESetType(pc, "boomeramg")
PetscOptionsSetValue("-pc_hypre_boomeramg_max_iter","2")
PetscOptionsSetValue("-pc_hypre_boomeramg_strong_threshold","0.25")
PetscOptionsSetValue("-pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi")
PCSetFromOptions(pc)
--]]

-- PROMETHEUS PC
---[[
PetscOptionsSetValue("-out_verbose","2") 
PCSetType(pc, "prometheus")
PCSetCoordinatesFromMesh(pc,mesh2)
PetscOptionsSetValue("-pc_mg_type","multiplicative")
PetscOptionsSetValue("-pc_mg_cycles","1")  
PetscOptionsSetValue("-prometheus_preduce_base","500") 
PetscOptionsSetValue("-prometheus_top_grid_limit","2500")
PetscOptionsSetValue("-aggmg_avoid_resmooth","")
PetscOptionsSetValue("-aggmg_smooths","1")
PetscOptionsSetValue("-prometheus_repartition","")
--PetscOptionsSetValue("-options_left","")
PCSetFromOptions(pc)
--]]
PetscOptionsPrint()

-- Solve system
KSPSetFromOptions(ksp)
--KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD)
KSPSolve(ksp,Fz,Uz)

--petsc_vector_dump(Uz,'tt2.txt')

rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)

-- Delete objects
KSPDestroy(ksp)
VecDestroy(Uz)
VecDestroy(Fz)
MatDestroy(Kz)
mesh2:delete()
mesh:delete()

print('Residual:',rnorm, ' NumIterations:', its)

te   = os.clock()
print(te-ts)
