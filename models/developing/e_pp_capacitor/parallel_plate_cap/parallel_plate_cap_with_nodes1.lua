-- Parallel plate capacitor test mesh

-- Uses Electrode element to tie all potentials on top and bottom electrode
--   % Top    potentials are short circuited to Node x1
--   % Bottom potentials are short circuited to Node x2

-- HiQLab
-- Copyright (c): Regents of the University of California

 
-- Set up properties

eps0  = 8.854e-13  -- Permittivity of free space
g     = 2e-6       -- Gap
w     = 50e-6      -- Gap width
dense = 1e-5       -- Mesh density
order = 1          -- Quadratic elements
meshtol = 1e-12;

mesh = Mesh:new(2);

-- External nodes to connect electrode
x1    = {25e-6, 7e-6}
x2    = {25e-6,-5e-6}
idx1  = mesh:add_node(x1);
idx2  = mesh:add_node(x2);
econ1 = {idx1}
econ2 = {idx2}

-- Parallel plate capacitor
D    = mesh:PMLScalar2d(eps0,0);
mesh:blocks2d({0,w}, {0,g}, D)

-- Add Electrode elements via global shape functions
idg1 = mesh:add_global(1)
idg2 = mesh:add_global(1)
ele1 = mesh:Electrode(idg1,1)
ele2 = mesh:Electrode(idg2,1)
eno1 = mesh:add_element(econ1, ele1, 1, 1)
eno2 = mesh:add_element(econ2, ele2, 1, 1)

-- Tie mesh together
mesh:tie()

-- Boundary conditions: bottom is at zero volts
mesh:set_bc(function(x,y)
  if y == g then return 'u', 0; end     -- Fix potentials that are 
  if y == 0 then return 'u', 0; end      --  connected to global variable

  if (x==x1[1]) and (y==x1[2]) 
            then return 'u', 1; end     -- Fix node x1 at potential of 1[V]
  if (x==x2[1]) and (y==x2[2]) 
            then return 'u', 0; end     -- Fix node x2 at potential of 0[V]
end)

-- Tie all the top potentials together with global shape function
function trial_shapeg1(x,y)
  if y == g then 
     return 1.0; 
  else 
     return 0.0;
  end
end

-- Tie all the bottom potentials together with global shape function
function trial_shapeg2(x,y)
  if y == 0 then 
     return 1.0; 
  else 
     return 0.0;
  end
end

-- Input for global shape functions must be a TABLE
mesh:set_globals{trial_shapeg1,trial_shapeg2}
