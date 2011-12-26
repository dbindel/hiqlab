
require 'common.lua'

-- Mesh parameters
order    = order    or 1                    -- Order of shape functions
dense    = dense    or order                -- Approximate node spacing
f0       = f0       or 0                    -- Stretch function parameters
bmLength = bmLength or 30e-6                -- Beam Length
bmWidth  = bmWidth  or  8e-6                -- Beam Width
stLength = stLength or  4e-6                -- Stub Length
stWidth  = stWidth  or  2e-6                -- Stub Width
wafer    = wafer    or  '100'               -- Wafer orientation
angle    = angle    or  0                   -- Crystal axis angle(rad)
-- etype    = etype    or  'planestrain'       -- Stress analysis type
dense    = dense*bmWidth
meshtol  = dense/1000
dense_x  = dense
dense_y  = dense
caxis    = { 0, 0, 1, 1, 0, 0}              -- Orientation of crystal axes

-- Circuit parameters
Vt_a = Vt_a  or 0.1  -- Top Plate
Vb_a = Vb_a  or   0  -- Bottom Plate

-- Nondimensionalize length parameters
beam_material = beam_material or 'aln'
dim_scales    = pz_nondim(beam_material,stWidth)
Vt_a,Vb_a     = material_normalize('V',Vt_a,Vb_a)

-- Construct mesh
mesh  = Mesh:new(2)
belt  = mesh:PMLElastic2d_pz_planestrain(beam_material, caxis)

x1 = -stLength - bmLength/2
x2 =           - bmLength/2
x3 = 0.0
x4 =-x2
x5 =-x1

y1 = -bmWidth/2
y2 = -stWidth/2
y3 =-y2
y4 =-y1


mesh:blocks2d( { x1, x2             }, 
               {     y2, y3         }, 
                belt, order, dense_x, dense_y )
mesh:blocks2d( {     x2, x3, x4     },
               { y1, y2, y3, y4     },
                belt, order, dense_x, dense_y )
--mesh:blocks2d( {             x4, x5 }, 
--               {     y2, y3         }, 
--                belt, order, dense_x, dense_y )
mesh:tie()

-- Define boundary condition

mesh:set_bc(function(x,y)
  if mesheq( x, x1) then return 'uu ', 0, 0; end
end)

function top_potential(x,y)
  if mesheq( y, y4) then return '  u',    1; end;
end;

function bot_potential(x,y)
  if mesheq( y, y1) then return '  u',    1; end;
end;

-- Tie some potentials together
ntie_t = mesh:TieField2(top_potential);
ntie_b = mesh:TieField2(bot_potential);
