require 'common.lua'

num_elem = num_elem or 6            -- Number of elements along block
order    = order    or 1            -- Element order
div_x    = order * num_elem
div_y    = order * num_elem * 6
meshtol  = 1e-6                     -- Mesh tolerance
Ri = 10                             -- Arch inner radius
Rt = 1                              -- Arch thickness

mesh  = Mesh:new(2);
etype = mesh:PMLElastic2d_planestrain{E = 1000.0, nu = 0.25}

-- Define mesh using block command
function ring(x,y)
   rnew = Ri + (x+1)/2 * Rt;
   anew =  0 + (y+1)/2 * pi; 
   return rnew*cos(anew),rnew*sin(anew);
end
mesh:add_block_transform(div_x+1, div_y+1, etype, order, ring)

-- Define boundary conditions
point_load(0, Ri+Rt, {uy = -5})
clamp_boundary(function(x,y) return mesheq(y,0) end, 'ux', 'uy')
