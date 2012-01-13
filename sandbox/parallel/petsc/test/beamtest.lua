require 'testcommon.lua'
b2d = loadfile('beam2d.lua')

mtype = mtype or 'seqbaij'
nbx   = nbx or 6*10
nby   = nby or 6*2

function save_results(text)
  return {text, r1, it1, r2, it2, r3, it3}
end

function print_saved(results)
  local text, r1, it1, r2, it2, r3, it3 = unpack(results)
  print('----' .. text .. '----')
  print_results('hi', 'hi', r1, it1)
  print_results('hi', 'lo', r2, it2)
  print_results('hi', 'lo', r3, it3)
end

------------------------------------------------------------------

experiments = {
 { '-- Second order, force loading --',
   { order = 2, force_load = true, nx = nbx/2, ny = nby/2 },
   { order = 1, force_load = true, nx = nbx,   ny = nby   } },
 { '-- Third order, force loading --',
   { order = 3, force_load = true, nx = nbx/3, ny = nby/3 },
   { order = 1, force_load = true, nx = nbx,   ny = nby   } },
 { '-- Second order, disp loading --',
   { order = 2, force_load = false,  nx = nbx/2, ny = nby/2 },
   { order = 1, force_load = false,  nx = nbx,   ny = nby   } },
 { '-- Third order, disp loading --',
   { order = 3, force_load = false,  nx = nbx/3, ny = nby/3 },
   { order = 1, force_load = false,  nx = nbx,   ny = nby   } },
 { '-- Second order, force loading (2x) --',
   { order = 2, force_load = true, nx = 2*nbx/2, ny = 2*nby/2 },
   { order = 1, force_load = true, nx = 2*nbx,   ny = 2*nby   } },
 { '-- Third order, force loading (2x) --',
   { order = 3, force_load = true, nx = 2*nbx/3, ny = 2*nby/3 },
   { order = 1, force_load = true, nx = 2*nbx,   ny = 2*nby   } },
 { '-- Second order, disp loading (2x) --',
   { order = 2, force_load = false,  nx = 2*nbx/2, ny = 2*nby/2 },
   { order = 1, force_load = false,  nx = 2*nbx,   ny = 2*nby   } },
 { '-- Third order, disp loading (2x) --',
   { order = 3, force_load = false,  nx = 2*nbx/3, ny = 2*nby/3 },
   { order = 1, force_load = false,  nx = 2*nbx,   ny = 2*nby   } } 
}

results = {}
for i,experiment in experiments do
  Khi, Fhi, Klo, Flo = form_problem(b2d, experiment[2], experiment[3])
  x1,r1,it1, x2,r2,it2, x3,r3,it3 = solve_problem(Khi, Fhi, Klo, Flo)
  clean_matrices(Khi, Klo)
  clean_vectors(x1, x2, x3, Fhi, Flo)
  results[i] = save_results(experiment[1])
end

for i in experiments do
  print_saved(results[i])
end
