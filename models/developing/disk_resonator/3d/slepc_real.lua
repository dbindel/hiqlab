ts   = os.clock()
-- Loading the Lua input file
f0    = 00
order = 1;
dense = 1e-6; -- 1e-6
--runfile = loadfile('../../lua_trilinos/test/block3d.lua')
runfile = loadfile('plate3d.lua')
--runfile = loadfile('pml1d.lua')

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()

mesh2 = mesh:get_clean_mesh()
mesh2:set_bc(bc_function)
mesh2:initialize()
mesh2:apply_bc()

numid = mesh2:get_numid()
print("Numid:",numid)
print("Tot:",mesh2:get_ndf()*mesh2:numnp())

-- Set KSP Options
PetscOptionsSetOption("-st_ksp_rtol","1.0e-5")
PetscOptionsSetOption("-st_ksp_type","gmres" )

-- Set PC Options
-- HYPRE PC
---[[
PetscOptionsSetOption("-st_pc_type", "hypre")
PetscOptionsSetOption("-st_pc_hypre_type", "boomeramg")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_max_iter","2")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_strong_threshold","0.25")
PetscOptionsSetOption("-st_pc_hypre_boomeramg_relax_type_all","symmetric-SOR/Jacobi")
PetscOptionsSetOption("-st_ksp_monitor","")
PetscOptionsPrint()
--]]
-- PROMETHEUS PC
--[[
PetscOptionsSetOption("-st_pc_type", "prometheus")
PetscOptionsSetOption("-st_ksp_monitor","")
--PetscOptionsSetValue("-log_summary","")
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
nev = 1
--w0  = 0
w0      = 2.889953e+08*2*pi;
w0      = 0;
status, dr, di, Vz = compute_eigs_slepc(mesh2, w0, nev)

-- Print eigenvalues
for i = 1,nev do
    local Q = sqrt(dr[i]^2+di[i]^2)/2/di[i]
    print(i, ':', dr[i]/2e6/pi, '  :', di[i]/2e6/pi, 'MHz',  '  Q:',Q)
end

-- Print results
--[[
ToMatrixMarketFile('Vr.mm',Vr)
ToMatrixMarketFile('Vi.mm',Vi)
--]]

-- Delete objects
VecDestroy(Vz)
mesh2:delete()

te   = os.clock()
print(te-ts)
