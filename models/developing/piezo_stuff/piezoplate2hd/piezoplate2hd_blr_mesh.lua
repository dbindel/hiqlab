
require 'common.lua'

-- Set default values for any unset parameters

order    = order    or 1                    -- Order of shape functions
dense    = dense    or order                -- Approximate node spacing
f0       = f0       or 20                   -- Stretch function parameters
bmLength = bmLength or  50e-6               -- Beam Length
bmWidth  = bmWidth  or 200e-6               -- Beam Width
pmlDepth = pmlDepth or  20e-6
anLength = anLength or  50e-6;
anWidth  = anWidth  or 100e-6; 
stLength = stLength or  95e-6                -- Stub Length
stWidth  = stWidth  or  8.0e-6                -- Stub Width
wafer    = wafer    or  '100'               -- Wafer orientation
angle    = angle    or  0                   -- Crystal axis angle(rad)
-- etype    = etype    or  '2hd'       -- Stress analysis type
dense    = dense*bmWidth
meshtol  = dense/1000
dense_x  = dense
dense_y  = dense
caxis    = { 1, 0, 0, 0, 1, 0}              -- Orientation of crystal axes

-- Nondimensionalize length parameters
beam_material = beam_material or 'aln'
dim_scales    = pz_nondim(beam_material,stWidth)

-- To pass values into matlab
beam_material_n = get_material(beam_material)
beam_material_n = pz_material_normalize(beam_material_n)
cL              = dim_scales.L
cT              = dim_scales.T
cM              = dim_scales.M
cA              = dim_scales.A

mesh  = Mesh:new(2)
belt  = mesh:PMLElastic2hd_pz(beam_material, caxis)

-- Mesh consists of two superblocks tied together
x5 = -stLength - bmLength/2 - anLength
x6 =  x5 + pmlDepth
x1 = -stLength - bmLength/2
x2 =           - bmLength/2
x3 = -x2
x4 = -x1
x7 = -x6
x8 = -x5

y1 = -anWidth/2;
y2 = -anWidth/2 + pmlDepth
y3 = -stWidth/2
y4 = -y3
y5 = -y2
y6 = -y1

mesh:blocks2d( {           x1,       x2    }, 
               {           y3,       y4    }, 
                belt, order, dense_x, dense_y )
mesh:blocks2d( {           x2,       x3    },
               { -bmWidth/2.0,       y3,       y4, bmWidth/2.0 },
                belt, order, dense_x, dense_y )
mesh:blocks2d( {           x3,       x4    }, 
               {           y3,       y4    }, 
                belt, order, dense_x, dense_y )
mesh:blocks2d( { x5, x6, x1 },
               { y1,y2,y3,y4,y5,y6},
                belt, order, dense_x, dense_y )
mesh:blocks2d( { x4, x7, x8 },
               { y1,y2,y3,y4,y5,y6},
                belt, order, dense_x, dense_y )
mesh:tie()

-- Define stretching function
belt:set_stretch(function(x,y)
  local xs = 0  -- x stretch value
  local ys = 0  -- y stretch value
  if meshleq(x,x6)  then xs = f0*(x6 - x )/(x6-x5) end
  if meshgeq(x,x7)  then xs = f0*(x  - x7)/(x8-x7) end
  if meshleq(y,y2)  and (meshleq(x,x1) or meshgeq(x,x4)) 
                    then ys = f0*(y2 - y )/(y2-y1) end
  if meshgeq(y,y5)  and (meshleq(x,x1) or meshgeq(x,x4))
                    then ys = f0*(y  - y5)/(y6-y5) end
  return xs, ys
end)

-- Define boundary condition

mesh:set_bc(function(x,y)
  if mesheq( x, x5) then return 'uuu', 0, 0, 0; end
  if mesheq(y, y1) and (meshleq(x,x1) or meshgeq(x,x4)) 
                    then return 'uuu', 0, 0, 0; end
  if mesheq( x, x8) then return 'uuu', 0, 0, 0; end
  if mesheq(y, y6) and (meshleq(x,x1) or meshgeq(x,x4)) 
                    then return 'uuu', 0, 0, 0; end
  if meshbetween(x, x2, x3) 
                    then return '  u', 1; end
  
  return '  u', 0;
              
end)

