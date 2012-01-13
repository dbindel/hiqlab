-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2);

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
-- etype = mesh:PMLElastic2d_planestrain(mtype);
etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define beam geometry
beam_l = 10
beam_w =  2

-- Divisions along edges
div_x  =  4
div_y  =  2

-- Define coordinates
x   = {}
for i = 0,div_x do
    for j = 0,div_y do
        x[2*(i*(div_y+1)+j)+1] = beam_l/div_x * i
        x[2*(i*(div_y+1)+j)+2] = beam_w/div_y * j
    end
end

-- Define element connectivity
con = {}
e_num = 0
for i = 0,div_x-1 do
    for j = 0,div_y-1 do
      e_num = e_num + 1
      con[4*e_num-3] = i   *(div_y+1) + j
      con[4*e_num-2] = i   *(div_y+1) + j+1
      con[4*e_num-1] =(i+1)*(div_y+1) + j 
      con[4*e_num  ] =(i+1)*(div_y+1) + j+1
    end
end

-- Add nodes and elements to mesh
mesh:add_node(x, (div_x+1)*(div_y+1))
mesh:add_element( con, etype, 4, e_num)

-- Define boundary conditions
meshtol = 1e-6
point_load(beam_l, 0, {uy = 1})
clamp_boundary(function(x,y) return mesheq(x,0,1e-6) end, 'ux', 'uy')