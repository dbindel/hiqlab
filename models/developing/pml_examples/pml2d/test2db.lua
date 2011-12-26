-- 2D test mesh

-- HiQLab
-- Copyright (c): Regents of the University of California


-- Set up properties

min = math.min
max = math.max
abs = math.abs

rho = 1;
E   = 10;
nu  = 0.3;

xrad = 12;  -- Radius of domain (x dimension)
yrad = 12;  -- Radius of domain (y dimension)
xpml =  9;  -- Start of PML (x dimension)
ypml =  9;  -- Start of PML (y dimension)
xsym =  0;  -- Symmetry boundary at x = 0?
mx   = 36;  -- Element spaces in x
my   = 24;  -- Element spaces in y
f0   = 40;  -- Maximum stretch parameter

w = 1;      -- Drive frequency
a = 1;      -- Magnitude of vertical pull in forcing
b = 0;      -- Magnitude of tilt in forcing
s = 1;      -- Stretch magnitude for displacement plot

order = 3;           -- Element order
nen   = (order+1)^2; -- Number of element nodes 

-- Assemble mesh quantities

mesh = Mesh:new(2);
D    = mesh:own( PMLElastic2d:new(E, nu, rho, 0) );

if xsym == 1 then
  mesh:add_block(0,-yrad, xrad,0, mx+1,my+1, D, order);
else
  mesh:add_block(-xrad,-yrad, xrad,0, mx+1,my+1, D, order);
end

D:set_stretch(function(x,y)
  local xs = max( (abs(x)-xpml)/(xrad-xpml) * f0/w, 0 )
  local ys = max( (abs(y)-ypml)/(yrad-ypml) * f0/w, 0 )
  return xs, ys
end)

