-- HiQLab
-- Copyright (c): Regents of the University of California


require 'common.lua'


-- Set default values for any unset parameters

order  = order  or  1              -- Order of shape functions
dense  = dense  or  order*0.5e-6   -- Approximate node spacing
f0     = f0     or 40              -- Stretch function parameter
rpost  = rpost  or 2.97e-6/2       -- Post radius
rpost2 = rpost2 or 3.22e-6/2       -- Post base radius
rdisk  = rdisk  or 40e-6           -- Disk radius
hdisk  = hdisk  or  1.5e-6         -- Height of disk layer
hpost  = hpost  or  1.05e-6        -- Height of post layer
rbd    = rbd    or  2e-6           -- Radius of domain
rpml   = rpml   or  4e-6           -- Radius of PML

post_material = post_material or 'sige'
disk_material = disk_material or 'sige'

meshtol = dense/100
mesh = Mesh:new(2)

pelt = make_material(post_material, 'axis')
delt = make_material(disk_material, 'axis')


-- Mesh consists of three superblocks tied together

function npts(l)  return 1+order*ceil( l/dense )  end

npostx = npts(rpost2)
nposty = npts(hpost)
ndiskx = npts(rdisk-rpost)
ndisky = npts(hdisk)

mesh:add_block_shape(ndiskx, ndisky, delt, order,
  {rpost,       hpost,   rpost,       hpost+hdisk, 
   rpost+rdisk, hpost,   rpost+rdisk, hpost+hdisk})
mesh:add_block_shape(npostx, ndisky, pelt, order, 
  {0,     hpost,   0,     hpost+hdisk, 
   rpost, hpost,   rpost, hpost+hdisk})
mesh:add_block_shape(npostx, ndisky, pelt, order, 
  {0,      0,   0,     hpost, 
   rpost2, 0,   rpost, hpost})
mesh:blocks2d( {0, rpost2, rbd, rpml},  {-rpml, -rbd, 0}, pelt, order, dense )

mesh:tie()

-- Define stretching function

pelt:set_stretch(function(x,y)
  local xs = 0   -- x stretch value
  local ys = 0   -- y stretch value
  if x >  rbd and y <= 0 then xs =  f0*(x-rbd)/(rpml-rbd) end
  if y < -rbd            then ys = -f0*(y+rbd)/(rpml-rbd) end
  return xs, ys
end)


-- Define boundary conditions

mesh:set_bc(function(x,y)
  if mesheq(x,0)                     then return 'u', 0; end
  if mesheq(x,rpost+rdisk) and y > 0 then return 'f', 1; end
end)
