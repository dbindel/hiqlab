require 'common.lua'

order   = order or 1  -- Order of element
dense   = 2.0e-6      -- Approximate element size
meshtol = dense/100   -- Default mesh tolerance

l  = 10e-6  -- Beam length
w  =  2e-6  -- Beam width
t  =  1e-6  -- Beam thickness
Pw = 20e-6  -- Anchor Width
Ph = 10e-6  -- Anchor Height
Pd =  2e-6  -- PML Depth
f0 = 0      -- Define PML paramter

-- Define element type
mesh = Mesh:new(2)
ted_nondim('silicon2',7e-6)
etype = mesh:PMLElastic2d_te_planestress('silicon2')

-- Define mesh and anchor using block commands
mesh:blocks2d( { 0, l}, { -w/2, w/2}, etype)
bc_func, st_func = pml_blocks2d({0,0},{w/2,-w/2},2,Pw,Ph,Pd,f0,
                                etype,order,dense,dense)
mesh:tie()

-- Set stretch function
etype:set_stretch(st_func)

-- Define boundary condition
function temperature_bc_function(x,y)
  if mesheq( y, w/2.0) and meshgeq(x, 0) then return '  u',        10; end
  if mesheq( y,-w/2.0) and meshgeq(x, 0) then return '  u',       -10; end
end
mesh:set_bc{bc_func,temperature_bc_function}

-- Define sensing function for tip displacement
tip_displacement2    = nodal2d_indicator('uy', l, -w/2)

