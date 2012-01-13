-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2);

-- Define element type
mtype     = {}
mtype.E   =  10
mtype.nu  =  0.3
mtype.rho =  1
etype = make_material_e(mesh, mtype, 'planestrain');

-- Define order of element
order = 3

-- Define mesh tolerance
meshtol = 1e-6

-- Define geometry of domain
xrad = 12  -- Radius of domain (x dimension)
yrad = 12  -- Radius of domain (y dimension)
xpml =  9  -- Start of PML (x dimension)
ypml =  9  -- Start of PML (y dimension)
mx   = 36  -- Element spaces in x
my   = 24  -- Element spaces in y

-- PML parameter
f0   = f0 or 40  -- Maximum stretch parameter

-- Forcing parameters
a = a or 1 -- Magnitude of vertical pull in forcing
b = b or 0 -- Magnitude of tilt in forcing

-- Define mesh using block command
mesh:add_block(-xrad,-yrad, xrad,0, mx+1,my+1, etype, order);

-- Define stretching function
function stretch_function(x,y)
  local xs = max( (abs(x)-xpml)/(xrad-xpml) * f0, 0 )
  local ys = max( (abs(y)-ypml)/(yrad-ypml) * f0, 0 )
  return xs, ys
end

etype:set_stretch(stretch_function)

-- Define boundary condition
function bc_function(x,y)
  if mesheq( y, 0) and meshleq( abs(x), 1) then return 'uu', 0, a+b*x; end
  if mesheq( y,-yrad)                      then return 'uu', 0, 0;     end
  if mesheq( abs(x), xrad)                 then return 'uu', 0, 0;     end
end

mesh:set_bc(bc_function)
