-- Include function definition file
require 'common.lua'

-- Construct and define dimension of mesh
mesh = Mesh:new(3)

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic3d(mtype);

-- Define block geometry
block_l = 10
block_w = 10
block_d = 10

-- Element order, approximate size of elements
order  =  3
dense  =  10

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks3d command
mesh:blocks3d({0, block_l},{0, block_w}, {0, block_d}, 
               etype, order, dense)

-- Define boundary conditions
meshtol = 1e-6

-- Define boundary conditions (displacement and force)
function bc_function(x,y,z)
  -- Displacement boundary conditions
  if mesheq( z,  0.0, 1e-6) then
                     return 'uuu', 0, 0, 0; 
  end
  -- Force boundary conditions
  if mesheq( z, 10.0, 1e-6) then
                     return '  f', 1.0; 
  end
end
mesh:set_bc(bc_function)

