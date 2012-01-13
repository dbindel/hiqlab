-- HiQLab
-- Copyright (c): Regents of the University of California
-- $Id: test_taxis.lua,v 1.2 2006/05/02 04:16:32 tkoyama Exp $

l = 10        -- Beam length
w = 2         -- Beam width
dense = 0.5   -- Approximate element size (for block generator)
order = 2     -- Order of elements
nen   = 9     -- Number of element nodes

E   = 1;      -- Young's modulus
nu  = 0.25;   -- Poisson ratio
rho = 1;      -- Mass density

ltype = ltype or -1;   -- Harmonic index in theta

if ltype < 0 then
  mesh  = Mesh:new(2)
  mat   = mesh:own( PMLElasticAxis:new(E, nu, rho) );
else
  mesh  = Mesh:new(2)
  mat   = mesh:own( PMLElasticTAxis:new(E, nu, rho, ltype) );
end

mesh:blocks2d( { 0, l }, { -w/2.0, w/2.0 }, mat, order, dense )

mesh:set_bc(function(x,y)
  if x == 0 then return 'uuu', 0, 0, 0; end
end)

