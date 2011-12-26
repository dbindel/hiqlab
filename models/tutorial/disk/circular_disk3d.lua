-- Include function definition file
require 'common.lua'
require 'circular_disk_func3d.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(3);

-- Define mesh tolerance
order    = order    or 1
dense   = 0.5      -- Approximate element size
meshtol = dense/100   -- Default mesh tolerance

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic3d(mtype);

-- Define geometric parameters
ri = 1;
rt1 = 5;
rt2 = rt1+dense*2;
rt3 = rt1+dense;
ro = 20;
rb1 = 3;
rb2 = 5;
rpml= 8;
nd = 10
xo = 1
yo = 0.0
t  = 1

-- Define number of elements along the block

-- Define mesh using block command
ang = 0
-- Disk portion
zt  = {4, 5}
mesh_inner_disk(xo,yo,ri,zt,nd,ang)
mesh_post_interface(xo,yo,ri,rt1,zt,nd,ang)
mesh_transition(rt1,rt3,rt2,zt,nd,ang)
mesh_steady(rt2,ro,zt,nd*3,ang)

-- Stem portion
zt  = {0, 4}
mesh_inner_disk(xo,yo,ri,zt,nd,ang)

-- Inside the block
zt  = {-rb2, 0}
mesh_inner_disk(xo,yo,ri,zt,nd,ang)
mesh_post_interface(xo,yo,ri,rb1,zt,nd,ang)
mesh_steady_end(rb1,rb2,zt,nd,ang)

-- Anchor
ndivpml = ceil((rpml-rb2)/dense)*order+1
ndivb   = ceil(rb2/dense)*order+1
zt  = {-rpml, -rb2}
mesh_inner_disk(xo,yo,ri,zt,nd,ang)
mesh_post_interface(xo,yo,ri,rb1,zt,nd,ang)
mesh_steady_end(rb1,rb2,zt,nd,ang)

mesh:blocks3dn( { rb2, rpml         }, 
                { ndivpml}, 
                {     0, rb2     },
                {nd+1}, 
                { -rpml, -rb2, 0 }, 
                { ndivpml, ndivb},
etype, order)
mesh:blocks3dn( 
                {     0, rb2     },
                {nd+1}, 
                { rb2, rpml         }, 
                { ndivpml}, 
                { -rpml, -rb2, 0 }, 
                { ndivpml, ndivb},
etype, order)
mesh:blocks3dn( 
                { rb2, rpml     },
                {ndivpml}, 
                { rb2, rpml         }, 
                { ndivpml}, 
                { -rpml, -rb2, 0 }, 
                { ndivpml, ndivb},
etype, order)

-- Tie mesh together
mesh:tie()

-- Define boundary conditions
--point_load(0,1, {uy = -5})
--clamp_boundary(function(x,y) return mesheq(x,0) end, 'ux')
--clamp_boundary(function(x,y) return mesheq(y,0) end, 'uy')
