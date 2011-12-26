-- 1D test mesh

-- HiQLab
-- Copyright (c): Regents of the University of California

 
-- Set up properties

min, max, abs = math.min, math.max, math.abs

xrad = 12;  -- Radius of domain (x dimension)
xpml =  9;  -- Start of PML (x dimension)
mx   = 36;  -- Element spaces in x
f0   = 20;  -- Maximum stretch parameter

nen   = order+1;  -- Number of element nodes 

-- Assemble mesh quantities

mesh = Mesh:new(1);
D    = mesh:own( PMLScalar1d:new(1, 1) );

mesh:add_block(0, xrad, mx+1, D, order);

D:set_stretch(function(x)
  return max( (abs(x)-xpml)/(xrad-xpml) * f0/k, 0 )
end)

-- Find where to apply BC

mesh:set_bc(function(x,y)
  if x == 0 then return 'u', 1  end
end)

