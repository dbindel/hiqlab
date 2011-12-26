-- Include function definition file
require 'common.lua'

-- Construct and define dimension of mesh
mesh = Mesh:new(2)

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define Coordinates
x  = { 
   0.0,  0.0,
   4.0,  0.0,
  10.0,  0.0,
   0.0,  4.5,
   5.5,  5.5,
  10.0,  5.0,
   0.0, 10.0,
   4.2, 10.0,
  10.0, 10.0 }

-- Element connectivity for linear quad elements
-- Explain how the connectivity is defined x-fix,y sweep
con = { 
   0, 3, 1, 4,
   1, 4, 2, 5,
   3, 6, 4, 7,
   4, 7, 5, 8}

-- Add nodes to mesh at once
-- Lua returns 0 based index for the node numbers
mesh:add_node(x,9);

-- Add elements to mesh
-- (connectivity, element type, number of element nodes, number of elements adding)
mesh:add_element(con, etype, 4, 4)

-- Define boundary conditions (displacement and force)
function bc_function(x,y)
  -- Displacement boundary conditions
  if mesheq( x,  0.0, 1e-6) and mesheq( y, 0.0, 1e-6) then
                     return 'uu', 0, 0; 
  end
  if mesheq( x,  0.0, 1e-6) then
                     return 'u ', 0;
  end
  -- Force boundary conditions
  if mesheq( x, 10.0, 1e-6) and mesheq( y, 5.0, 1e-6) then
                     return 'f ', 5.0; 
  end
  if mesheq( x, 10.0, 1e-6) then
                     return 'f ', 2.5;
  end
end
mesh:set_bc(bc_function)

