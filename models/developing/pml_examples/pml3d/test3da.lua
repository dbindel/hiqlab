-- 3D test mesh

-- HiQLab
-- Copyright (c): Regents of the University of California

 
-- Set up properties

min, max, abs = math.min, math.max, math.abs

rho = 1;
E   = 10;
nu  = 0.3;

xrad = 10;  -- Radius of domain (x dimension)
zrad = 10;  -- Radius of domain (z dimension)
xpml =  6;  -- Start of PML (x dimension)
zpml =  6;  -- Start of PML (z dimension)
mx   = 10;  -- Element spaces in x
mz   =  4;  -- Element spaces in z
f0   = 20;  -- Maximum stretch parameter

w = 1;      -- Drive frequency
a = 1;      -- Magnitude of vertical pull in forcing
b = 0;      -- Magnitude of tilt in forcing
s = 1;      -- Stretch magnitude for displacement plot

order = 2;           -- Element order
nen   = (order+1)^3; -- Number of element nodes 

-- Assemble mesh quantities

mesh = Mesh:new(3);
D    = mesh:own( PMLElastic3d:new(E, nu, rho) );

mesh:add_block(-xrad,-xrad,-xrad, xrad,xrad,0, mx+1,mx+1,mz+1, D, order);

D:set_stretch(function(x,y,z)
  local xs = max( (abs(x)-xpml)/(xrad-xpml) * f0/w, 0 )
  local ys = max( (abs(y)-xpml)/(xrad-xpml) * f0/w, 0 )
  local zs = max( (abs(z)-zpml)/(zrad-zpml) * f0/w, 0 )
  return xs, ys, zs
end)

-- Find where to apply BC

if nil then
mesh:set_bc(function(x,y,z)
  if z == 0 and abs(x) < 1e-8 and abs(y) < 1e-8 then 
    return 'uu', 0, a+b*x  
  end
end)
end

