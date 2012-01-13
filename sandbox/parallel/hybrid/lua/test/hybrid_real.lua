ts   = os.clock()
-- Loading the Lua input file
f0 = 00
order = 2;
dense = 1e-6;
--dense = 5e-6;
runfile = loadfile('plate3d.lua')
--runfile = loadfile('../../lua_trilinos/test/block3d.lua')
--runfile = loadfile('pml1d.lua')
--runfile = loadfile('disk3d.lua')

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

-- Assemble distributed stiffness(Petsc)
mattype    = 3;
is_reduced = 2;
scaleK       = 1e2
scaleM       = 1e16
--w0         = 2.889953e+08*2*pi;
--w0         = 1.5891068
--w0         = 2.2085182261522
--w0      = 2.789953e+10*2*pi/sqrt(scaleM*scaleK)
--w0      = 2.989953e+1*2*pi
w0           = 0
K = mesh2:assemble_dR_petsc(1, 0, -w0*w0*scaleM/scaleK, mattype, is_reduced)
MatScale(K,scaleK)

-- Assemble distributed mass(trilinos)
M = mesh2:assemble_dR_trilinos(0,0,1,1)
M:Scale(scaleM)

-- Construct KSP for stiffness
ksp = KSPCreate()
KSPSetOperators(ksp,K,K,"DIFFERENT_NONZERO_PATTERN")
KSPSetTolerances(ksp, 1.0e-9, PETSC_DEFAULT, PETSC_DEFAULT, PETSC_DEFAULT)

-- Set KSP Options
PetscOptionsSetOption("-ksp_rtol","1.0e-9")
PetscOptionsSetOption("-ksp_type","gmres" )
PetscOptionsSetOption("-ksp_monitor","")
--PetscOptionsSetValue("-log_summary","")
PetscOptionsPrint()

-- Set PC Options
-- HYPRE PC
--[[
PetscOptionsSetOption("-pc_type", "hypre")
PetscOptionsSetOption("-pc_hypre_type", "boomeramg")
PetscOptionsSetOption("-pc_hypre_boomeramg_max_iter","2")
PetscOptionsSetOption("-pc_hypre_boomeramg_strong_threshold","0.25")
PetscOptionsSetOption("-pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi")
PetscOptionsPrint()
--]]
-- PROMETHEUS PC
---[[
PetscOptionsSetValue("-out_verbose","2") 
PetscOptionsSetOption("-pc_type", "prometheus")
PetscOptionsSetOption("-pc_mg_type","multiplicative")
PetscOptionsSetOption("-pc_mg_cycles","2")  
PetscOptionsSetOption("-prometheus_preduce_base","500") 
PetscOptionsSetOption("-prometheus_top_grid_limit","2500")
PetscOptionsSetOption("-aggmg_avoid_resmooth","")
PetscOptionsSetOption("-aggmg_smooths","1")
PetscOptionsSetOption("-prometheus_repartition","")
PetscOptionsSetOption("-mg_levels_ksp_type","fgmres")
PetscOptionsSetOption("-mg_levels_pc_type","gs")
--PetscOptionsSetOption("-options_left","")
PetscOptionsPrint()
--]]

-- Construct Preconditioning operator to pass to Anasazi
pop = Petsc_Operator:new(ksp,M)
pop:SetFromOptions()
pop:SetMesh(mesh2)
--KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD)

-- Set Anasazi parameters
PL = ParameterList:new()
PL:set_int   ("Solver", Anasazi_BlockKrylovSchur)
PL:set_int   ("Block Size",        2)
PL:set_int   ("Max Blocks",        4)
PL:set_int   ("Max Restarts",    300)
PL:set_double("Convergence Tolerance",          1.0e-9)
--PL:set_double("Tol",          1.0e-9)
PL:set_string("Sigma",          "LM")    -- Sort manager parameters
PL:set_int("Verbosity", Anasazi_Warning+Anasazi_FinalSummary+Anasazi_IterationDetails+Anasazi_Debug+Anasazi_OrthoDetails) -- Output manager parameters

-- Set memory for eigenvector
PL:set_int("IsSymmetric", 0)
nev= 1
Vr = Epetra_MultiVector:new(numid,1*nev)

-- Compute eigenvalues
dr = {}
di = {}

status = compute_eigs_anasazi1(
pop,
nev,
PL, dr, di, Vr)
-- Print eigenvalues
for i = 1,nev do
   print(dr[i],' ',di[i])
end
undo_spectral_trans(nev,w0,0,1,dr,di)

-- Print eigenvalues
print('\n')
for i = 1,nev do
   print(dr[i],' ',di[i])
   print(dr[i]/2/pi,' ',di[i]/2/pi)
end

-- Print results
--[[
ToMatrixMarketFile('Vr.mm',Vr)
ToMatrixMarketFile('Vi.mm',Vi)
--]]


-- Delete objects
Vr:delete()
PL:delete()
pop:delete()
KSPDestroy(ksp)
M:delete()
MatDestroy(K)
mesh2:delete()
mesh:delete()

te   = os.clock()
print(te-ts)
