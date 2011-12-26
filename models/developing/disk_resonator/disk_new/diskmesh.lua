-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define mesh related global parameters
order  = order or  1              -- Order of shape functions
dense  = dense or  order*0.5e-6   -- Approximate node spacing
meshtol= dense / 100

-- Define element type
post_material = post_material or 'silicon2'
disk_material = disk_material or 'silicon2'
pelt = mesh:PMLElasticAxis(post_material)
delt = mesh:PMLElasticAxis(disk_material)

-- Set default values for any unset parameters
f0    = f0    or 40              -- Stretch function parameter
rpost = rpost or  1e-6           -- Post radius
rdisk = rdisk or 10e-6           -- Disk radius
hdisk = hdisk or  2e-6           -- Height of disk layer
hpost = hpost or  5e-7           -- Height of post layer
rbd   = rbd   or  2e-6           -- Radius of domain
rpml  = rpml  or  4e-6           -- Radius of PML

-- Mesh consists of three superblocks tied together
mesh:blocks2d( {0, rpost, rbd, rpml},  {-rpml, -rbd, 0}, pelt, order, dense )
mesh:blocks2d( {0, rpost}, {0, hpost, hpost + hdisk},    pelt, order, dense )
mesh:blocks2d( {rpost, rdisk},  {hpost, hpost + hdisk},  pelt, order, dense )

-- Tie mesh together
mesh:tie()

-- Define stretching function
function stretch_function(x,y)
  local xs = 0   -- x stretch value
  local ys = 0   -- y stretch value
  if x > rbd and meshleq(y,0) then xs =  f0*(x-rbd)/(rpml-rbd) end
  if y <-rbd                  then ys = -f0*(y+rbd)/(rpml-rbd) end
  return xs, ys
end
pelt:set_stretch(stretch_function)
--delt:set_stretch(stretch_function)

-- Define boundary conditions
function bc_function(x,y)
  if mesheq(x,0)                            then return 'u', 0; end
  if mesheq(x,rdisk) and       y > 0        then return 'f', 1; end
end
mesh:set_bc(bc_function)

-- Define driving and sensing functions
function bode_drive_function(x,y)
  if mesheq(x,rdisk) and  y > 0        then return 'f', 1; end
end
function bode_sense_function(x,y)
  if mesheq(x,rdisk) and  y > 0        then return 'u', 1; end
end

-- Indicator versions
bode_drive_func2 = make_nodal_func(bode_drive_function, 'f')
bode_sense_func2 = make_nodal_func(bode_sense_function, 'u')
