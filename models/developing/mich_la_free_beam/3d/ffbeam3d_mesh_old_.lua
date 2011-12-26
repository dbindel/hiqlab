
require 'common.lua'


-- Set default values for any unset parameters

-- Order and density of elements
order    = order    or 1               -- Order of shape functions
dense    = dense    or order*1.0       -- Approximate node spacing


-- Beam Parameters
bmLength = bmLength or 20e-6           -- Beam length
bmWidth  = bmWidth  or 2.0e-6          -- Beam width
bmHeight = bmHeight or 2e-6            -- Beam height
tabLength=tabLength or 9e-6
ftLength = ftLength or 25.6e-6      
ftWidth  = ftWidth  or 1.2e-6

-- Substrate Parameters
subWidth  = subWidth  or 20e-6  -- Width of substrate
subLength = subLength or 20e-6  -- Height of substrate
subDepth  = subDepth  or 10e-6  -- Depth of substrate
pmlDepth  = pmlDepth  or  4e-6  -- Depth of PML
f0        = f0        or 20     -- Stretch function parameter

-- Anchor Parameters
anWidth  =  anWidth   or  5e-6  -- Width of Anchor
anLength =  anLength  or  5e-6  -- Length of Anchor
anHeight =  anHeight  or  4e-6  -- Height of Anchor

-- Wafer Orientation
wafer    = wafer    or '100'           -- Wafer orientation
angle    = angle    or 0               -- Crystal axis angle(radians)

-- Element type
etype_b  = etype_b   or '3d'           -- Stress analysis type of beam
etype_s  = etype_s   or '3d'           -- Stress analysis type of substrate

-- Readjust density according to bmWidth
dense    = dense*bmWidth
meshtol  = dense/1000
dense_x  = dense
dense_y  = dense
dense_z  = dense

-- Nondimensionalize length parameters
beam_material  = beam_material or 'poly_silicon'
sub_material   = sub_material  or 'poly_silicon'
dim_scales     = ted_nondim(beam_material,bmWidth)
--beam_material = beam_material or 'sc_silicon'
--sub_material  = sub_material  or 'sc_silicon'

-- Passing values to matlab
beam_material_n=get_material(beam_material)
beam_material_n=ted_material_normalize(beam_material_n)
T0             = beam_material_n.T0
cL             = dim_scales.L
cT             = dim_scales.T
cM             = dim_scales.M
cTh            = dim_scales.Th


-- Dimensions for mesh

blRatio  = 0.25;                        -- Boundary Layer to Beam Center 
                                        --  Portion ratio
bcLength = bmLength*(1 - 2*blRatio);    -- Beam Center length
blLength = bmLength*blRatio;            -- Boundary Layer length
--anWidth  = bmWidth * 6.0;               -- Anchor Width
--anDepth  = bmWidth * 3.0;               -- Anchor Length
--ftWidth  = 0.60 * bmWidth;              -- Side footing length
--ftLength = 12*bmWidth;                  -- Side footing height
--tabLength= 0.5*bmLength                 -- Beam tab length
--pmlDepth = anDepth / 2.0;               -- PML Depth

-- Coordinates for mesh
x1 = - tabLength - ftWidth - bmLength/2.0;
x2 =             - ftWidth - bmLength/2.0;
x3 =                       - bmLength/2.0;
x4 = 0.0;

xs1= -bmLength/2.0 - ftWidth/2.0  - subWidth/2.0;
xs2= xs1 + pmlDepth;
xs3= -bmLength/2.0 - ftWidth/2.0 - anWidth/2.0;
xs4= x2;
xs5= x3;
xs6= -bmLength/2.0 - ftWidth/2.0 + anWidth/2.0;
xs8= -bmLength/2.0 - ftWidth/2.0 + subWidth/2.0;
xs7= xs8 - pmlDepth;

y1 = - bmWidth/2.0 - ftLength;
y2 = - bmWidth/2.0;
y3 = 0.0;
y4 = -y2;
y5 = -y1;

ys1= - bmWidth/2.0 - ftLength - anLength/2.0 - subLength/2.0;
ys2= ys1 + pmlDepth;
ys3= - bmWidth/2.0 - ftLength - anLength;
ys4= y1;
ys6= - bmWidth/2.0 - ftLength - anLength/2.0 + subLength/2.0;
ys5= ys6 - pmlDepth;

ys7 = -ys6;
ys8 = -ys5;
ys9 = -ys4;
ys10= -ys3;
ys11= -ys2;
ys12= -ys1; 

z1 = -subDepth;
z2 = -subDepth + pmlDepth;
z3 = 0.0;
z4 = anHeight - bmHeight;
z5 = anHeight;

-- Generate mesh

mesh  = Mesh:new(3)
elt_b = make_material_te(beam_material, etype_b)
elt_s = make_material_te( sub_material, etype_s)
--elt_b = make_material_te(beam_material, etype_b, wafer, angle)
--elt_s = make_material_te( sub_material, etype_s, wafer, angle)

-- Mesh consists of superblocks tied together

-- Main beam
mesh:blocks3d( { x1, x2         }, {     y2, y3 ,y4     },{ z4, z5 }, elt_b, 
                                 order, dense/3.0, dense, dense )--1
mesh:blocks3d( {         x3, x4 }, {     y2, y3 ,y4     },{ z4, z5 }, elt_b,
                                 order, dense/3.0, dense, dense )--2
mesh:blocks3d( {     x2, x3     }, {     y2, y3 ,y4     },{ z4, z5 }, elt_b,
                                 order, dense    , dense, dense )--3

-- Foot beams(Density 3 times)
mesh:blocks3d( {     x2, x3     }, { y1, y2             },{ z4, z5 }, elt_b,
                                 order, dense    , dense/3.0, dense )--4
mesh:blocks3d( {     x2, x3     }, {             y4, y5 },{ z4, z5 }, elt_b,
                                 order, dense    , dense/3.0, dense )--5

-- Anchor
mesh:blocks3d( { xs3, xs4, xs5, xs6 }, { ys3, ys4  }, { z3, z4, z5}, elt_b,
                                 order, dense    , dense,     dense )--6
mesh:blocks3d( { xs3, xs4, xs5, xs6 }, { ys9, ys10 }, { z3, z4, z5}, elt_b,
                                 order, dense    , dense,     dense )--8

-- Substrate
mesh:blocks3d( { xs1, xs2, xs3, xs4 , xs5 , xs6, xs7, xs8 },
          { ys1, ys2, ys3, ys4 , ys5 , ys6  },
          {  z1,  z2,  z3 }, elt_s,
                                 order, dense    , dense,     dense )--7
mesh:blocks3d( { xs1, xs2, xs3, xs4 , xs5 , xs6, xs7, xs8 },
          { ys7, ys8, ys9, ys10, ys11, ys12 },
          {  z1,  z2,  z3 }, elt_s,
                                 order, dense    , dense,     dense )--9

mesh:tie()

-- Define stretching function
elt_s:set_stretch(function(x,y,z)
  local xs = 0  -- x stretch value
  local ys = 0  -- y stretch value
  local zs = 0  -- z stretch value

  local x0 = - bmLength/2.0 - ftWidth/2.0;
  local y0m= - bmWidth/2.0  - ftLength - anLength/2.0;
  local y0p=   bmWidth/2.0  + ftLength + anLength/2.0;
  local z0 = 0.0;
  local noPmlx = subWidth/2.0  - pmlDepth;
  local noPmly = subLength/2.0 - pmlDepth;
  local noPmlz = subDepth      - pmlDepth;
 
  if meshleq(z,0) then
    -- Pml in the x direction
    if meshgeq(math.abs(x - x0),noPmlx) then
        xs = f0*(math.abs(x-x0 )-noPmlx) / pmlDepth
    end
    -- Pml in the y direction
    if meshleq(y,y2)  then
      if meshgeq(math.abs(y - y0m),noPmly) then
        ys = f0*(math.abs(y-y0m)-noPmly) / pmlDepth
      end
    elseif meshgeq(y,y4) then
      if meshgeq(math.abs(y - y0p),noPmly) then
        ys = f0*(math.abs(y-y0p)-noPmly) / pmlDepth
      end
    end
    -- Pml in the z direction
    if meshgeq(math.abs(z - z0),noPmlz) then
        zs = f0*(math.abs(z-z0 )-noPmlz) / pmlDepth
    end
  end
  return xs, ys, zs
end)

-- Define boundary conditions
-- Both disp and thermal bc
mesh:set_bc(function(x,y,z)

  -- Symmetric Boundary condition
  if meshbetween(y,y2,y4) and mesheq(x,0.0) and meshbetween(z,z4,z5) then
      return 'u   ', 0;
  end;

  if meshleq(z,0) then
    if meshleq(y,y2) or meshgeq(y,y4) then
        -- Fixed Edge Boundary conditions
        if mesheq(x,xs1) or mesheq(x,xs8) or mesheq(y,ys1) 
        or mesheq(y,ys6) or mesheq(y,ys7) or mesheq(y,ys12)
                         or mesheq(z,z1)  then
            return 'uuuu', 0, 0, 0, 0;
        end;
        -- Thermal boundary conditions in the interior
        if meshleq(x,xs2)  or  meshgeq(x,xs7)   then return '   u',      0; end
        if meshleq(y,ys2)  or  meshgeq(y,ys11)  then return '   u',      0; end
        if meshbetween(y,ys5,ys8)               then return '   u',      0; end
        if meshleq(z,z2)                        then return '   u',      0; end
    end
  end
end)
