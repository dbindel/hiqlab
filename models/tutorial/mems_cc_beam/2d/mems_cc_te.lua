-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
ted_nondim('sc_silicon',2e-6)

-- Define mesh related global parameters
order    = order    or 1               -- Order of shape functions
dense    = dense    or order*1e-6      -- Approximate node spacing
dense_x  = dense
dense_y  = dense
meshtol  = dense/100

-- Define geometry of domain
bmL = bmL or 14e-6           -- Beam length
bmW = bmW or 2e-6            -- Beam width
fcW = fcW or 4e-6            -- Forcing Width

-- Define element type
wafer          = wafer    or '100'           -- Wafer orientation
angle          = angle    or 0               -- Crystal axis angle(radians)
-- etype          = etype    or 'planestrain'   -- Stress analysis type
beam_material  = beam_material or 'sc_silicon'
belt  = mesh:PMLElastic2d_te_planestrain(beam_material, wafer, angle)

-- Mesh consists of one superblock tied together
mesh:blocks2d( { -bmL/2.0, 0, bmL/2.0 }, 
               { -bmW/2.0   , bmW/2.0  }, belt, order, dense_x, dense_y )

-- Tie mesh together
mesh:tie()

-- Define boundary conditions
function bc_function(x,y)
  if mesheq(y, -bmW/2.0) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
  if mesheq(y,  bmW/2.0) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
  
  if mesheq(x, -bmL/2.0)  then return 'uuu', 0, 0, 0; end
  if mesheq(x,  bmL/2.0)  then return 'uuu', 0, 0, 0; end
end

mesh:set_bc(bc_function)
