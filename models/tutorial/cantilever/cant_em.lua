require 'common.lua'

order   = order or 1       -- Order of element
dense   = dense or 2.0e-6  -- Approximate element size
meshtol = dense/100        -- Default mesh tolerance

l  = 10e-6     -- Beam length
w  =  2e-6     -- Beam width
t  =  1e-6     -- Beam thickness
g  =  2e-6     -- Gap width
P  =  3e-6     -- Force at tip
V  = V or 100  -- Voltage applied

mesh = Mesh:new(2)
em_nondim('silicon2',2e-6)
etype_e  = mesh:PMLElastic2d_planestrain('silicon2')
etype_em = mesh:CoupleEM2d(dim_scales.eps)

-- Define mesh geometry
mesh:blocks2d( { 0, l}, {  g, g+w}, etype_e)
mesh:blocks2d( { 0, l}, {  0,   g}, etype_em, order, dense, g)
mesh:tie()

-- Define boundary condition
clamp_boundary(
  function(x,y) return mesheq(x,0) or y < g end, 
  'ux', 'uy')

function gap_v_bc(x,y)
  if mesheq(y, 0) then return '  u',  0; end
  if mesheq(y, g) then return '  u',  V; end
end

mesh:set_bc{gap_v_bc}

-- Define drive and sense functions
tip_displacement2    = nodal2d_indicator('uy', l, g)
bode_sense_function2 = nodal2d_indicator('uy', l, g)
bode_force_function2 = nodal2d_indicator('uy', l, g, P/t)
