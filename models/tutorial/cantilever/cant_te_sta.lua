require 'common.lua'

order   = order or 1  -- Order of element
dense   = 2.0e-6      -- Approximate element size
meshtol = dense/100   -- Default mesh tolerance

l  = 10e-6  -- Beam length
w  =  2e-6  -- Beam width
t  =  1e-6  -- Beam thickness

-- Define element type
mesh = Mesh:new(2)
ted_nondim('silicon2',7e-6)
etype = mesh:PMLElastic2d_te_planestress('silicon2')

-- Define mesh geometry
mesh:blocks2d( { 0, l}, { -w/2, w/2}, etype)
mesh:tie()

-- Define boundary condition
function bc_function(x,y)
  if mesheq(y, w/2) and mesheq(x, 0) then return 'uuu', 0, 0,  10; end
  if mesheq(y,-w/2) and mesheq(x, 0) then return 'uuu', 0, 0, -10; end
  if mesheq(y, w/2) then return '  u',        10; end
  if mesheq(y,-w/2) then return '  u',       -10; end
  if mesheq(x,   0) then return 'uu ', 0, 0;      end
end

mesh:set_bc(bc_function)

-- Define sensing function for tip displacement
function tip_displacement(x,y)
  if mesheq(x, l) and mesheq(y, -w/2) then return ' u ', 1; end
end 

tip_displacement2 = make_nodal_func(tip_displacement, 'u')
