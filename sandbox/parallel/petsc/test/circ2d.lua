require 'common.lua'

E     = E or 1000
nu    = nu or 0.25
rho   = rho or 1
l     = l or 10
w     = w or 2
order = order or 3
nr    = nr or 10
nt    = nt or 10

mesh = Mesh:new(2);
etype = mesh:PMLElastic2d(E, nu, rho, 0);
-- mesh:add_block(0, -w/2, l, w/2, order*nx+1, order*ny+1, etype, order);
mesh:add_block_transform(order*nr+1, order*nt+1, etype, order,
  function(xx,yy)
    local theta = yy*pi/2
    local r     = (xx+2)/2
    return r*cos(theta), r*sin(theta)
  end);

-- Define boundary conditions
mesh:set_bc(function(x,y)
  if mesheq( x, 0, 1e-6) then return 'uu', 0, 0; end
  if mesheq( y, 0, 1e-6) and mesheq(x, 1.5, 1e-6) then return 'u ', 1; end
end)

