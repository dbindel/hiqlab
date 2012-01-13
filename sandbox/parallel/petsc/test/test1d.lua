-- 1D test mesh

-- HiQLab
-- Copyright (c): Regents of the University of California
 
-- Set up properties

order = order or 3    -- Polynomial order
nelt  = nelt or 1000  -- Number of elements

-- Assemble mesh quantities

mesh = Mesh:new(1);
D1   = mesh:PMLScalar1d(1, 1);
mesh:add_block(0, 1, nelt*order+1, D1, order);

-- Define BC

mesh:set_bc(function(x)  
  if x == 0          then return 'u', 1 end  
  if abs(x-1) < 1e-8 then return 'u', 0 end
end)

