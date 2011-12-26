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
belt  = mesh:PMLELastic2d_te_planestrain(beam_material, wafer, angle)

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
anW = bmW *12.0;          -- Anchor Width
anH = bmW * 6.0;          -- Anchor Length
pmlD= anH/ 2.0;           -- PML Depth

-- Coordinates for mesh
x1 = - bmL/2.0
x2 = - bcL/2.0
x3 = -x2
x4 = -x1
y1 = - bmW/2.0
y2 = -y1

-- Mesh consists of superblock tied together
mesh:blocks2d( { x1, x2, x3, x4}, { y1, y2}, belt, order, dense_x, dense_y )

bc_func1,st_func1 = pml_blocks2d({x1,x1},{y2,y1},3,anW,anH,pmlD,f0,
                                  belt,order,dense_x,dense_y)
bc_func2,st_func2 = pml_blocks2d({x4,x4},{y1,y2},3,anW,anH,pmlD,f0,
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
  if mesheq(y, -bmW/2.0) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
  if mesheq(y,  bmW/2.0) and meshbetween(x,-fcW/2.0,fcW/2.0) then
                                    return ' f ',    1   ; end
end
