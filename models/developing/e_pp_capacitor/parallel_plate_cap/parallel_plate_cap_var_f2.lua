-- Parallel plate capacitor test mesh(variant)

-- Uses TieField2 Element to tie 
--   1. Bottom    portion of potentials
--   2. Top left  portion of potentials 
--   3. Top right portion of potentials

-- HiQLab
-- Copyright (c): Regents of the University of California

 
-- Set up properties

eps0  = 8.854e-13  -- Permittivity of free space
g     = 20e-6       -- Gap
w     = 50e-6      -- Gap width
dense = 1e-6       -- Mesh density
order = 1          -- Quadratic elements
meshtol = 1e-12;

-- Assemble mesh quantities

x0    =  0
x1    = 20e-6
x2    = 30e-6
x3    = w

mesh = Mesh:new(2);
D    = mesh:PMLScalar2d(eps0,0);
mesh:blocks2d({x0,x1,x2,x3}, {0,g}, D)
mesh:tie()

-- Boundary conditions

mesh:set_bc(function(x,y)
--  Must set blank boundary conditions
end)

-- Tie all the bottom potentials together
ntie1 = mesh:TieField2(function(x,y)
  if mesheq(y,0) then return 'u', 1; end
end)
-- Tie all the top left  potentials together
ntie2 = mesh:TieField2(function(x,y)
  if mesheq(y,g) and meshleq(x,w/2-1e-6) then return 'u', 1; end
end)
-- Tie all the top right potentials together
ntie3 = mesh:TieField2(function(x,y)
  if mesheq(y,g) and meshgeq(x,w/2+1e-6) then return 'u', 1; end
end)
