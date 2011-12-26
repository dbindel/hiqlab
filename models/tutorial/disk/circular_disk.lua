-- Include function definition file
require 'common.lua'
require 'circular_disk_func.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2);

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define mesh tolerance
dense   = 0.5      -- Approximate element size
meshtol = dense/100   -- Default mesh tolerance
meshtol = 1.0e-6

-- Define geometric parameters
ri = 1;
rt1 = 5;
rt2 = rt1+dense*2;
rt3 = rt1+dense;
ro = 20;
nd = 2
xo = 1
yo = 0.0

-- Define number of elements along the block
num_elem = num_elem or 1
order    = order    or 1
div_x    = order*num_elem
div_y    = order*num_elem
tran     = 1

-- Define mesh using block command
ang = 0
mesh_inner_disk(xo,yo,ri,nd,ang)
mesh_post_interface(xo,yo,ri,rt1,nd,ang)
mesh_transition(rt1,rt3,rt2,nd,ang)
mesh_steady(rt2,ro,nd*3,ang)
--mesh_steady_end(rt2,ro,nd*3,ang)

-- Tie mesh together
mesh:tie()

-- Define boundary conditions
point_load(0,1, {uy = -5})
clamp_boundary(function(x,y) return mesheq(x,0) end, 'ux')
clamp_boundary(function(x,y) return mesheq(y,0) end, 'uy')
