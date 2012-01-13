-- Default parameters --

mtype     = mtype     or 'seqaij'     -- PETSc matrix type
nelt      = nelt      or 1000         -- Number of elements in mesh
hi_order  = hi_order  or 3            -- Order of higher-order mesh

-----------------------------------------------------
require '../petsc-util.lua'
t1d = loadfile('test1d.lua')

print "-- Assemble high-order problem --"
order = hi_order;
t1d()
mesh:initialize()
K3 = mesh:assemble_dR_petsc(1,0,0, mtype)
F3 = mesh:assemble_R_petsc()
mesh:delete()

print "-- Assemble linear problem ---"
nelt = hi_order*nelt;
order = 1;
t1d()
mesh:initialize()
K1 = mesh:assemble_dR_petsc(1,0,0, mtype)
F1 = mesh:assemble_R_petsc()

print "-- Test solver --"
opt = {mesh=use_coord and mesh, rtol=1e-6, maxits=5000}
x1,r1,it1 = petsc_solve_x{K3, F3, opt}
x2,r2,it2 = petsc_solve_x{K3, F3, opt, B=K1}
x3,r3,it3 = petsc_solve_x{K1, F1, opt}

print "-- Destroy matrices --"
petsc_destroy_vec(x1)
petsc_destroy_vec(x2)
petsc_destroy_vec(x3)
petsc_destroy_mat(K1)
petsc_destroy_mat(K3)
mesh:delete()

print "-- Done --"
function print_results(hi_order, lo_order, r, i)
  print("(" .. hi_order .. ", " .. lo_order .. "): ", r, i)
end
print_results(hi_order, hi_order, r1, it1)
print_results(hi_order, 1,        r2, it2)
print_results(1,        1,        r3, it3)
