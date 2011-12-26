-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
ted_nondim('sc_silicon',2e-6)

-- Define mesh related global parameters
order    = order    or 1               -- Order of shape functions
dense    = dense    or order*1e-6      -- Approximate node spacing
dense_x  = dense
dense_y  = dense
meshtol  = dense/100

-- Define element type
wafer          = wafer    or '100'           -- Wafer orientation
angle          = angle    or 0               -- Crystal axis angle(radians)
-- etype          = etype    or 'planestrain'   -- Stress analysis type
beam_material  = beam_material or 'sc_silicon'
belt  = mesh:PMLElastic2d_te_planestrain(beam_material, wafer, angle)

-- Set PML parameters
f0       = f0       or 20               -- Stretch function parameter

-- Define geometry of domain
bmL = bmL or 14e-6        -- Beam length
bmW = bmW or 2e-6         -- Beam width
fcW = fcW or 4e-6         -- Forcing Width

blR = 0.25;               -- Boundary Layer to Beam Center
                                        --  Portion ratio
bcL = bmL*(1 - 2*blR);    -- Beam Center length
blL = bmL*blR;            -- Boundary Layer length
anW = bmW * 6.0;          -- Anchor Width
anH = bmW * 3.0;          -- Anchor Length
pmlD= anH/ 2.0;           -- PML Depth

ftW  = 0.60 * bmW         -- Side footing length
ftH  = 12   * bmW         -- Side footing height
tabL = 0.5  * bmL         -- Beam tab length

-- Coordinates for mesh
x1 = - tabL - ftW - bmL/2.0;
x2 =        - ftW - bmL/2.0;
x3 =              - bmL/2.0;
x4 = -x3;
x5 = -x2;
x6 = -x1;

y1 = 0.0;
y2 = ftH;
y3 = ftH + bmW;

-- Mesh consists of 3 superblocks tied together

mesh:blocks2d( { x1, x2, x3, x4, x5, x6 }, { y2, y3 }, belt, order, dense_x, dense_y )
mesh:blocks2d( {     x2, x3             }, { y1, y2 }, belt, order, dense_x, dense_y )
mesh:blocks2d( {             x4, x5     }, { y1, y2 }, belt, order, dense_x, dense_y )

bc_func1,st_func1 = pml_blocks2d({x2,x3},{y1,y1},3,anW,anH,pmlD,f0,
                                  belt,order,dense_x,dense_y)
bc_func2,st_func2 = pml_blocks2d({x4,x5},{y1,y1},3,anW,anH,pmlD,f0,
                                  belt,order,dense_x,dense_y)

-- Tie mesh together
mesh:tie()

-- Define stretching function
belt:set_stretch{st_func1,st_func2}

-- Define boundary conditions
-- Both disp and no thermal bc
mesh:set_bc{bc_func1,bc_func2}

-- Sensing function
function bode_force_pattern(x,y)
  -- Forced boundary conditions
  if mesheq(y, y2) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
  if mesheq(y, y3) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
end
