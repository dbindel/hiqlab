-----------------------------------------------------
-- 
-- Khi, Fhi, Klo, Flo = load_problem(load_name)
-- Khi, Fhi, Klo, Flo = form_problem(meshfile, hi_env, lo_env)
-- x1,r1,it1, x2,r2,it2, x3,r3,it3 = solve_problem(Khi, Fhi, Klo, Flo)
-- save_problem(Khi, Fhi, Klo, Flo, save_name)
-- save_results(x1, x2, x3, save_name)
-- print_results(hi_order, lo_order, r, i)
-- clean_vectors(...)
-- clean_matrices(...)
-- 
-----------------------------------------------------

require '../petsc-util.lua'

rtol = rtol or 1e-6;        -- relative tolerance for solves


function set_env(e)
  local saved_env = {}
  for key,val in e do
    saved_env[key] = _G[key]
    _G[key] = val
  end
  return saved_env
end


function form_problem(meshfile, hi_env, lo_env, keep_mesh)

  -- Assemble higher-order problem --
  local saved_env = set_env(hi_env)
  meshfile()
  mesh:initialize()
  local Khi = mesh:assemble_dR_petsc(1,0,0, mtype)
  local Fhi = mesh:assemble_R_petsc()
  mesh:delete()
  set_env(saved_env)

  -- Assemble linear problem ---
  saved_env = set_env(lo_env)
  meshfile()
  mesh:initialize()
  local Klo = mesh:assemble_dR_petsc(1,0,0, mtype)
  local Flo = mesh:assemble_R_petsc()
  if not keep_mesh then mesh:delete() end
  set_env(saved_env)

  return Khi, Fhi, Klo, Flo, mesh

end


function load_problem(load_name)

  -- Load problem --
  local Khi = petsc_matrix_load(load_name .. 'Khi')
  local Fhi = petsc_vector_load(load_name .. 'Fhi')
  local Klo = petsc_matrix_load(load_name .. 'Klo')
  local Flo = petsc_vector_load(load_name .. 'Flo')

  return Khi, Fhi, Klo, Flo

end


function solve_problem(Khi, Fhi, Klo, Flo, mesh)

  -- Test solver --
  local opt = { rtol=rtol, mesh=mesh, maxits=5000 }
  local x1,r1,it1 = petsc_solve_x{Khi, Fhi, opt}
  local x2,r2,it2 = petsc_solve_x{Khi, Fhi, opt, B=Klo}
  local x3,r3,it3 = petsc_solve_x{Klo, Flo, opt}

  return x1,r1,it1, x2,r2,it2, x3,r3,it3

end


function save_problem(Khi, Fhi, Klo, Flo, save_name)

  -- Save matrices --
  petsc_matrix_bdump(Khi, save_name .. 'Khi')
  petsc_vector_bdump(Fhi, save_name .. 'Fhi')
  petsc_matrix_bdump(Klo, save_name .. 'Klo')
  petsc_vector_bdump(Flo, save_name .. 'Flo')

end


function save_results(x1, x2, x3, save_name)

  -- Save result vectors --
  petsc_vector_save(x1,  save_name .. 'x1')
  petsc_vector_save(x2,  save_name .. 'x2')
  petsc_vector_save(x3,  save_name .. 'x3')

end


function clean_vectors(...)
  for i = 1,arg.n do
    petsc_destroy_vec(arg[i])
  end
end


function clean_matrices(...)
  for i = 1,arg.n do
    petsc_destroy_mat(arg[i])
  end
end


function print_results(hi_order, lo_order, r, i)
  print("(" .. hi_order .. ", " .. lo_order .. "): ", r, i)
end
