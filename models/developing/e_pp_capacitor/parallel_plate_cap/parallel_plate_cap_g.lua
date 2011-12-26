-- Parallel plate capacitor test mesh

-- Uses global shape functions to tie all potentials on top electrode

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
  if y == g then return 'u', 0; end -- For global shape function
                                    -- Must make boundary value 0
end)

-- Tie all the top potentials together with global shape function
idg1 = mesh:add_global(1)
function trial_shapeg1(x,y)
  if y == g then 
     return 1.0; 
  else 
     return 0.0;
  end
end

-- Input for global shape functions must be a TABLE
mesh:set_globals{trial_shapeg1}
