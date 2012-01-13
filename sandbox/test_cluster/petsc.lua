-- Start clock
ts   = os.clock()

-- Loading the Lua input file
runfile = loadfile('testg.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()

-- Get clean mesh(for Prometheus)
---[[
mesh2 = mesh:get_clean_mesh()
mesh2:set_bc(bc_function)
mesh2:initialize()
mesh2:apply_bc()
--]]

numid  = mesh2:get_numid()
mypid  = MPI_Comm_rank(PETSC_COMM_WORLD)
numproc= MPI_Comm_size(PETSC_COMM_WORLD)
if (mypid==0) then
  print("Number of processes:",numproc)
  print("NDOF:",numid)
  -- Open file to write data
  filename = string.format("petsc_results_%d.txt",numproc)
  f = io.open(filename,"w")
  f:write('Number of processes:',numproc,'\n')
  f:write("NDOF:",numid,'\n')
end

-- Assemble distributed stiffness
mattype    = 3; -- 2(MPIAIJ), 3(MPIBAIJ)
is_reduced = 2; -- 2(IDENTITY_DIRICHLET), 1(REMOVE_DIRICHLET)
Kz = mesh2:assemble_dR_petsc(1, 0, 0, mattype, is_reduced)

-- Assemble the RHS
Fz = mesh2:assemble_R_petsc(is_reduced)
VecScale(Fz,-1)

-- Assemble LHS
Uz = MatGetVecX(Kz)
VecSetFromOptions(Uz)

-- Create KSP
ksp = KSPCreate()
KSPSetOperators(ksp, Kz, Kz, "DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-9, PETSC_DEFAULT, PETSC_DEFAULT, 2000)
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
--PetscOptionsSetValue("-out_verbose","2") 
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

-- Solve system
--PetscOptionsPrint()
KSPSetFromOptions(ksp)
ts_s   = os.clock()
KSPSolve(ksp,Fz,Uz)
te_s   = os.clock()

-- Check residual
rnorm = KSPGetResidualNorm(ksp)
its   = KSPGetIterationNumber(ksp)
R     = MatGetVecY(Kz)
MatMult(Kz, Uz, R)
VecAXPY(R, -1, Fz)
rnorm2 = qVecNorm(R,2)
rnormf = qVecNorm(Fz,2)

-- Write file for plotting
--[[
print('MyPID:',mypid)
Mesh_SetU_Petsc_Vec(mesh2,Uz,is_reduced)
if mypid==0 then
  dxf = DXFile:new('testg0')
  dxf:writemesh(mesh2)
  dxf:delete()
end
--]]

-- Delete objects
KSPDestroy(ksp)
VecDestroy(Uz)
VecDestroy(Fz)
MatDestroy(Kz)
mesh2:delete()
mesh:delete()

if (mypid==0) then
  print('Norm of residual(prec)  :',rnorm)
  print('Norm of residual(abs)   :',rnorm2)
  print('Norm of residual(rel)   :',rnorm2/rnormf)
  print('Number of Iterations    :',its)
  f:write('Norm of residual(prec)  :',rnorm,'\n')
  f:write('Norm of residual(abs)   :',rnorm2,'\n')
  f:write('Norm of residual(rel)   :',rnorm2/rnormf,'\n')
  f:write('Number of Iterations    :',its,'\n')
end

-- Stop clock
te   = os.clock()
if (mypid==0) then
   print('Elapsed time            :',te-ts)
   print('Elapsed time(Only solve):',te_s-ts_s)
   f:write('Elapsed time            :',te-ts,'\n')
   f:write('Elapsed time(Only solve):',te_s-ts_s,'\n')
   f:close()
end
