-- HiQLab
-- Copyright (c): Regents of the University of California

-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define mesh related global parameters
order  = order  or  1              -- Order of shape functions
dense  = dense  or  order*0.5e-6   -- Approximate node spacing
meshtol= dense/100

-- Define element type
post_material = post_material or 'sige'
disk_material = disk_material or 'sige'
pelt = mesh:PMLElasticAxis(post_material)
delt = mesh:PMLELasticAxis(disk_material)

-- Set default values for any unset parameters
f0     = f0     or 40              -- Stretch function parameter
rpost  = rpost  or 2.97e-6/2       -- Post radius
rpost2 = rpost2 or 3.22e-6/2       -- Post base radius
rdisk  = rdisk  or 40e-6           -- Disk radius
hdisk  = hdisk  or  1.5e-6         -- Height of disk layer
hpost  = hpost  or  1.05e-6        -- Height of post layer
rbd    = rbd    or  2e-6           -- Radius of domain
rpml   = rpml   or  4e-6           -- Radius of PML


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
mesh:blocks2d( { 0, rpost2, rbd, rpml },  { -rpml, -rbd, 0 }, pelt )

-- Tie mesh together
mesh:tie()

-- Define stretching function
function stretch_function(x,y)
  local xs = 0   -- x stretch value
  local ys = 0   -- y stretch value
  if meshsg(x, rbd) and meshleq(y,0) then xs =  f0*(x-rbd)/(rpml-rbd) end
  if meshsl(y,-rbd)                  then ys = -f0*(y+rbd)/(rpml-rbd) end
  return xs, ys
end
pelt:set_stretch(stretch_function)

-- Define boundary conditions
function bc_function(x,y)
  if mesheq(x,0)                            then return 'u', 0; end
  if mesheq(x,rpost+rdisk) and meshsg(y,0)  then return 'f', 1; end
end
mesh:set_bc(bc_function)

-- Define forcing function
function force_pattern(x,y)
  if mesheq(x,rpost+rdisk) and meshsg(y,0)  then return 'f', 1; end
end
