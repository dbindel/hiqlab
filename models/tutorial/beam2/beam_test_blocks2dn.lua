-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2);

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define beam geometry
beam_l = 10
beam_w =  2

-- Element order, divisions along edges
order  =  1
div_x  =  4
div_y  =  2

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks2dn command
mesh:blocks2dn({0, beam_l},{div_x+1},{0, beam_w},{div_y+1}, etype, order)

-- Define boundary conditions
meshtol = 1e-6
point_load(beam_l, 0, {uy = 1})
clamp_boundary(function(x,y) return mesheq(x,0,1e-6) end, 'ux', 'uy')
