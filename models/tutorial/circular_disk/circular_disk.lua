-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2);

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define number of elements along the block
num_elem = num_elem or 1
order    = order    or 1
div_x    = order*num_elem
div_y    = order*num_elem

-- Define mesh tolerance
meshtol = 1e-6

-- Define nodal coordinates
xc = {
  0.0, 0.0,
  0.0, 0.5,
  0.5, 0.0,
  0.4, 0.4
}
xr1= {
  0.5, 0.0,
  0.4, 0.4,
  1.0, 0.0,
  cos(pi/4),sin(pi/4)
}
curv1 = { 0, 1, 0, 0}
xr2= {
  0.4, 0.4,
  0.0, 0.5,
  cos(pi/4),sin(pi/4),
  0.0, 1.0
}
curv2 = { 0, 1, 0, 0}

-- Define mesh using block command
mesh:add_block_shape(div_x+1, div_y+1, etype, order, xc)
mesh:add_curved_block_shape2d(div_x+1, div_y+1, etype, order, xr1, curv1)
mesh:add_curved_block_shape2d(div_x+1, div_y+1, etype, order, xr2, curv2)

-- Tie mesh together
mesh:tie()

-- Define boundary conditions
point_load(0,1, {uy = -5})
clamp_boundary(function(x,y) return mesheq(x,0) end, 'ux')
clamp_boundary(function(x,y) return mesheq(y,0) end, 'uy')
