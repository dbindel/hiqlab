require 'common.lua'

order   = 3     -- Define order of element
meshtol = 1e-6  -- Define mesh tolerance

xrad = 12  -- Radius of domain (x dimension)
yrad = 12  -- Radius of domain (y dimension)
xpml =  9  -- Start of PML (x dimension)
ypml =  9  -- Start of PML (y dimension)
mx   = 36  -- Element spaces in x
my   = 24  -- Element spaces in y

f0 = f0 or 40  -- Maximum stretch parameter

a = a or 1     -- Magnitude of vertical pull in forcing
b = b or 0     -- Magnitude of tilt in forcing

mesh = Mesh:new(2);
etype = mesh:PMLElastic2d_planestrain{ E = 10, nu = 0.3, rho = 1 }

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
