-- Include function definition file
require 'common.lua'

-- Construct and define dimension of mesh
mesh = Mesh:new(3)

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
mtype.rho =   10.0
scale     =   1e-6
--mech_nondim(mtype,scale)
etype = make_material_e( mtype, '3d');

-- Define block geometry
block_l = 10*scale
block_w = 10*scale
block_d = 10*scale

-- Element order, approximate size of elements
order  =  1
dense  =  10.0*scale

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks3d command
mesh:blocks3d({0, block_l},{0, block_w}, {0, block_d}, 
               etype, order, dense)

-- Define boundary conditions
meshtol = scale*1e-4

-- Define boundary conditions (displacement and force)
function bc_function(x,y,z)
  -- Displacement boundary conditions
  if mesheq( z,  0.0, meshtol) then
                     return 'uuu', 0, 0, 0; 
  end
  -- Force boundary conditions
  if mesheq( z, block_d, meshtol) then
                     return '  f', 1.0; 
  end
end
mesh:set_bc(bc_function)

