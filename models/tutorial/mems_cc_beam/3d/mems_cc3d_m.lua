-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(3)

-- Define mesh related global parameters
order    = order    or 1               -- Order of shape functions
dense    = dense    or order*1e-6      -- Approximate node spacing
meshtol  = dense/100

-- Define geometry of domain
bmL = bmL or 14e-6           -- Beam length
bmW = bmW or 2e-6            -- Beam width
bmT = bmT or 2e-6            -- Beam thickness
fcW = fcW or 4e-6            -- Forcing Width

-- Define element type
wafer          = wafer    or '100'           -- Wafer orientation
angle          = angle    or 0               -- Crystal axis angle(radians)
-- etype          = '3d'                        -- Stress analysis type
beam_material  = beam_material or 'sc_silicon'
belt  = mesh:PMLElastic3d_te(beam_material, wafer, angle)

-- Mesh consists of one superblocks tied together
mesh:blocks3d( { -bmL/2.0, 0, bmL/2.0 }, 
               { -bmW/2.0   , bmW/2.0  }, 
               { -bmT/2.0   , bmT/2.0 }, 
               belt, order, dense)

-- Tie mesh together
mesh:tie()

-- Define boundary conditions
function bc_function(x,y,z)
  if mesheq(x, -bmL/2.0)  then return 'uuu', 0, 0, 0; end
  if mesheq(x,  bmL/2.0)  then return 'uuu', 0, 0, 0; end
end

mesh:set_bc(bc_function)
