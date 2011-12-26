
require 'common.lua'

-- Set default values for any unset parameters

order    = order    or 1                    -- Order of shape functions
dense    = dense    or order                -- Approximate node spacing
f0       = f0       or 20                   -- Stretch function parameters
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

-- Anchor paramters
pmlDepth = pmlDepth or   5e-6
anLength = anLength or  15e-6;
anWidth  = anWidth  or  30e-6;

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
belt  = mesh:PMLElastic2d_pz_planestrain(beam_material, caxis)

x1 = -stLength - bmLength/2.0
x2 =           - bmLength/2.0
x3 = -x2
x4 = -x1
y1 = -bmWidth/2.0
y2 = -stWidth/2.0
y3 = -y2
y4 = -y1

xa1 = x1 - anLength
xa2 = xa1 + pmlDepth
xa3 = x1
xa4 = -xa3
xa5 = -xa2
xa6 = -xa1

ya1 = -anWidth/2
ya2 = ya1 + pmlDepth
ya3 = y2
ya4 =-ya3
ya5 =-ya2
ya6 =-ya1

-- Mesh consists of two superblocks tied together

mesh:blocks2d( {          x1, x2        }, 
               {          y2, y3        }, 
                belt, order, dense_x, dense_y )
mesh:blocks2d( {              x2, x3    },
               {      y1, y2, y3, y4    },
                belt, order, dense_x, dense_y )
mesh:blocks2d( {                  x3, x4}, 
               {          y2, y3        }, 
                belt, order, dense_x, dense_y )
mesh:blocks2d( { xa1, xa2, xa3 },
               { ya1, ya2, ya3, ya4, ya5, ya6},
                belt, order, dense_x, dense_y )
mesh:blocks2d( { xa4, xa5, xa6},
               { ya1, ya2, ya3, ya4, ya5, ya6},
                belt, order, dense_x, dense_y )

mesh:tie()

-- Define stretching function
belt:set_stretch(function(x,y)
  local xs = 0  -- x stretch value
  local ys = 0  -- y stretch value
  if meshleq(x,xa2)  then xs = f0*(xa2 - x )/pmlDepth end
  if meshgeq(x,xa5)  then xs = f0*(x  - xa5)/pmlDepth end
  if meshleq(y,ya2)  and (meshleq(x,xa3) or meshgeq(x,xa4))
                    then ys = f0*(ya2 - y )/pmlDepth end
  if meshgeq(y,ya5)  and (meshleq(x,xa3) or meshgeq(x,xa4))
                    then ys = f0*(y  - ya5)/pmlDepth end
  return xs, ys
end)

-- Define boundary condition

mesh:set_bc(function(x,y)
  if mesheq(x, xa1) then return 'uu ',  0, 0; end
  if mesheq(y, ya1) and (meshleq(x,xa3) or meshgeq(x,xa4))
                    then return 'uu ',  0, 0; end
  if mesheq(x, xa6) then return 'uu ',  0, 0; end
  if mesheq(y, ya6) and (meshleq(x,xa3) or meshgeq(x,xa4))
                    then return 'uu ',  0, 0; end
  if mesheq(y, y1) and meshbetween(x, x2, x3)
                    then return '  u', 1; end
  if mesheq(y, y4) and meshbetween(x, x2, x3)
                    then return '  u',-1; end

end)

