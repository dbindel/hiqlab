require 'common.lua'

E     = E or 1000
nu    = nu or 0.25
rho   = rho or 1
l     = l or 10
w     = w or 2
order = order or 3
nx    = nx or 50
ny    = ny or 10

mesh = Mesh:new(2);
etype = mesh:PMLElastic2d(E, nu, rho, 0);
mesh:add_block(0, -w/2, l, w/2, order*nx+1, order*ny+1, etype, order);

-- Define boundary conditions
mesh:set_bc(function(x,y)
  if mesheq( x, 0, 1e-6) then return 'uu', 0, 0; end
  if mesheq( x, l, 1e-6) and mesheq( y, 0, 1e-6) then 
    if force_load then return ' f', 1 else return ' u', 1 end
  end
end)

