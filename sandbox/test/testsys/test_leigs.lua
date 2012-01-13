beamf = loadfile 'beammesh.lua'

-- Load mesh and compute eigs
beamf()
mesh:initialize()
sigma = 1e16
verbose = 0
K = mesh:assemble_dR(1, 0, 0)
M = mesh:assemble_dR(0, 0, 1)
Ks1 = mesh:assemble_dR(1, 0, -sigma)
Ks2 = mesh:assemble_dR(1, 0,  sigma)

function test_eig(d, dknown, msg)
  assert(abs((d(1)-dknown)/dknown) < 1e-4,
         msg .. ': Incorrect eigenvalue')
end


---------------------------
-- Test symmetric driver --
---------------------------

-- Mode 1: Standard driver
d,v = leigs{
  op = K,
  nev = 1,
  symmetric = true,
  mode = 'standard',
  verbose = verbose
}
test_eig(d, 1.4694e+12, 'Symmetric standard problem')

-- Mode 2: Generalized
d,v = leigs{
  op = function(x,y)
    K:apply(x,y)
    x:copy(y)
    M:solve(y,x) 
  end,
  M = M,
  nev = 1,
  real = true,
  symmetric = true,
  which = 'SM',
  mode = 'generalized',
  maxitr = 300,
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Symmetric generalized problem (standard)')

-- Mode 3: Shift-invert
d,v = leigs{
  op = K,
  M = M,
  nev = 1,
  symmetric = true,
  mode = 'shift-invert',
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Symmetric generalized problem (shift-invert)')

-- Mode 4: Buckling
d,v = leigs{
  op = function(x,y)
    Ks1:solve(y,K:apply(x))
  end,
  M = K,
  sigma = sigma,
  nev = 1,
  symmetric = true,
  mode = 'buckling',
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Symmetric generalized problem (buckling)')

-- Mode 5: Cayley
d,v = leigs{
  op = function(x,y)
    Ks1:solve(y,Ks2:apply(x))
  end,
  M = M,
  sigma = sigma,
  nev = 1,
  symmetric = true,
  mode = 'Cayley',
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Symmetric generalized problem (Cayley)')

------------------------------
-- Test nonsymmetric driver --
------------------------------

-- Mode 1: Standard driver
d,v = leigs{
  op = K,
  nev = 1,
  mode = 'standard',
  verbose = verbose
}
test_eig(d, 1.4694e+12, 'Unsymmetric standard problem')

-- Mode 2: Generalized
d,v = leigs{
  op = function(x,y)
    K:apply(x,y)
    x:copy(y)
    M:solve(y,x) 
  end,
  M = M,
  nev = 1,
  real = true,
  which = 'SM',
  mode = 'generalized',
  maxitr = 300,
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Unsymmetric generalized problem (standard)')

-- Mode 3: Shift-invert
d,v = leigs{
  op = K,
  M = M,
  nev = 1,
  mode = 'shift-invert',
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Unsymmetric generalized problem (shift-invert)')


------------------------------
-- Test complex driver      --
------------------------------

-- Mode 1: Standard driver
d,v = leigs{
  op = K,
  nev = 1,
  real = false,
  mode = 'standard',
  verbose = verbose
}
test_eig(d, 1.4694e+12, 'Complex standard problem')

-- Mode 2: Generalized
d,v = leigs{
  op = function(x,y)
    K:apply(x,y)
    x:copy(y)
    M:solve(y,x) 
  end,
  M = M,
  nev = 1,
  real = false,
  which = 'SM',
  mode = 'generalized',
  maxitr = 300,
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Complex generalized problem (standard)')

-- Mode 3: Shift-invert
d,v = leigs{
  op = K,
  M = M,
  nev = 1,
  real = false,
  mode = 'shift-invert',
  verbose = verbose
}
test_eig(d, 3.0748e+16, 'Complex generalized problem (shift-invert)')
