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
l =  2
w =  2

-- Element order, approximate size of elements
order  =  1
dense_x=  2
dense_y=  2
meshtol= 1e-6

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks2d command
mesh:blocks2d({0, l},{ 0, w}, etype, order, dense_x, dense_y)
mesh:blocks2d({0, l},{-w, 0}, etype, order, dense_x, dense_y)

-- Tie mesh together
mesh:tie()

function bc_function(x,y)
    if mesheq(y,-w) then return 'uu', 0, 0; end;
    if mesheq(y, w) then return 'f ', 1   ; end;
end
mesh:set_bc(bc_function)

-- Define boundary conditions
--point_load(beam_l, 0, {uy = 1})
--clamp_boundary(function(x,y) return mesheq(x,0,1e-6) end, 'ux', 'uy')
