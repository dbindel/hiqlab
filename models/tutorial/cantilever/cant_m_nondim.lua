require 'common.lua'

order   = order or 1  -- Order of element
dense   = 2.0e-6      -- Approximate element size
meshtol = dense/100   -- Default mesh tolerance

l  = 10e-6  -- Beam length
w  =  2e-6  -- Beam width
t  =  1e-6  -- Beam thickness
P  =  3e-6  -- Define force at tip

-- Define element type
mesh = Mesh:new(2)
mech_nondim('silicon2',7e-6)
etype = mesh:PMLElastic2d_planestress('silicon2')

-- Define mesh geometry
mesh:blocks2d( { 0, l}, { -w/2, w/2}, etype)
mesh:tie()

-- Define boundary condition
function bc_function(x,y)
  if mesheq(x, 0) then return 'uu', 0, 0; end
  if mesheq(x, l) and mesheq(y, -w/2) then return ' f',  P/t; end
end

mesh:set_bc(bc_function)

-- Define sensing function for tip displacement
tip_displacement2    = nodal2d_indicator('uy', l, -w/2)
bode_sense_function2 = nodal2d_indicator('uy', l, -w/2)
bode_force_function2 = nodal2d_indicator('uy', l, -w/2, P/t)
