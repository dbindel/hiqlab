ts   = os.clock()
-- Loading the Lua input file
f0 = 00
order = 2;
dense = 1;
runfile = loadfile('plate3d.lua')
--runfile = loadfile('../../lua_trilinos/test/block3d.lua')
--runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print("Numid:",numid)

-- Set KSP Options
PetscOptionsSetOption("-st_ksp_rtol","1.0e-9")
PetscOptionsSetOption("-st_ksp_type","gmres" )
PetscOptionsSetOption("-st_ksp_monitor","")
--PetscOptionsSetValue("-log_summary","")

-- Set PC Options
-- HYPRE PC
---[[
PetscOptionsSetOption("-st_pc_type", "hypre")
PetscOptionsSetOption("-st_pc_hypre_type", "boomeramg")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_max_iter","2")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_strong_threshold","0.25")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi")
PetscOptionsPrint()
--]]
-- PROMETHEUS PC
--[[
PetscOptionsSetOption("-st_pc_type", "prometheus")
--PetscOptionsSetValue("-out_verbose","2") 
PetscOptionsSetOption("-pc_mg_type","multiplicative")
PetscOptionsSetOption("-pc_mg_cycles","1")  
PetscOptionsSetOption("-prometheus_preduce_base","500") 
PetscOptionsSetOption("-prometheus_top_grid_limit","2500")
PetscOptionsSetOption("-aggmg_avoid_resmooth","")
PetscOptionsSetOption("-aggmg_smooths","1")
PetscOptionsSetOption("-prometheus_repartition","")
--PetscOptionsSetOption("-options_left","")
PetscOptionsPrint()
--]]








-- Assemble distributed stiffness
nev = 2
w0  = 0
status, dr, di, Vz = compute_eigs_slepc(mesh, w0, nev)

-- Print eigenvalues
for i = 1,nev do
   print(dr[i],' ',di[i])
end

-- Print results
--[[
ToMatrixMarketFile('Vr.mm',Vr)
ToMatrixMarketFile('Vi.mm',Vi)
--]]

-- Delete objects
VecDestroy(Vz)
mesh:delete()

te   = os.clock()
print(te-ts)
