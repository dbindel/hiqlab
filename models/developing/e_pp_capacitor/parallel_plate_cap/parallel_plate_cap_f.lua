-- Parallel plate capacitor test mesh

-- Uses TieField Element to tie all potentials on top electrode

-- HiQLab
-- Copyright (c): Regents of the University of California

 
-- Set up properties

eps0  = 8.854e-13  -- Permittivity of free space
g     = 2e-6       -- Gap
w     = 50e-6      -- Gap width
dense = 1e-5       -- Mesh density
order = 1          -- Quadratic elements
meshtol = 1e-12;

-- Assemble mesh quantities

mesh = Mesh:new(2);
D    = mesh:PMLScalar2d(eps0,0);
mesh:blocks2d({0,w}, {0,g}, D)
mesh:tie()

-- Boundary conditions: bottom is at zero volts

mesh:set_bc(function(x,y)
  if y == 0 then return 'u', 0; end
end)

-- Tie all the top potentials together

ntie = mesh:TieField(function(x,y)
  if mesheq(y,g) then return 'u', 1; end
end)
