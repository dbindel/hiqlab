-- Include function definition file
require 'common.lua'

-- Construct and define dimension of mesh
mesh = Mesh:new(3)

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = make_material_e( mtype, '3d');
mtype1    = {}
mtype1.E   =  500.0
mtype1.nu  =    0.25
etype1 = make_material_e( mtype1, '3d');

-- Define block geometry
block_x = 10
block_y = 10
block_z = 10
plate_t =  1
arm_l   = 10
arm_w   =  2

-- Element order, approximate size of elements
order  =  1
dense  =  0.15

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks3d command
mesh:blocks3d({-block_x/2,-block_x/2+arm_w,  block_x/2-arm_w, block_x/2},
              {-block_y/2,-block_y/2+arm_w,  block_y/2-arm_w, block_y/2}, 
			  {0, block_z-plate_t, block_z}, 
               etype, order, dense)
mesh:blocks3d({ block_x/2, block_x/2+arm_l},
              { block_y/2-arm_w, block_y/2},
			  { block_z-plate_t, block_z  },
			  etype1, order, dense)
mesh:blocks3d({-block_x/2,-block_x/2+arm_w},
              { block_y/2, block_y/2+arm_l},
			  { block_z-plate_t, block_z  },
			  etype1, order, dense)
mesh:blocks3d({-block_x/2-arm_l,-block_x/2},
              {-block_y/2,-block_y/2+arm_w},
			  { block_z-plate_t, block_z  },
			  etype1, order, dense)
mesh:blocks3d({ block_x/2-arm_w, block_x/2},
              {-block_y/2-arm_l,-block_y/2},
			  { block_z-plate_t, block_z  },
			  etype1, order, dense)

-- Define boundary conditions
meshtol = 1e-6
mesh:tie()

-- Define boundary conditions (displacement and force)
function bc_function(x,y,z)
  -- Displacement boundary conditions
  if mesheq( z,  0.0, 1e-6) then
                     return 'uuu', 0, 0, 0; 
  end
  -- Force boundary conditions
  if mesheq( abs(x), block_x/2+arm_l, 1e-1) then
                              return '  f', 1.0; 
  end
  if mesheq( abs(y), block_y/2+arm_l, 1e-1) then
                              return '  f', 1.0; 
  end
end
mesh:set_bc(bc_function)

