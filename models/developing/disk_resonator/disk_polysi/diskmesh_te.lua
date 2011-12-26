-- HiQLab
-- Copyright (c): Regents of the University of California

-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
mtype = 'sc_silicon'
--mtype = 'poly_silicon'
ted_nondim(mtype,2e-6)

-- Define mesh related global parameters
order  = order or  1              -- Order of shape functions
dense  = dense or  order*0.5e-6   -- Approximate node spacing
meshtol= dense / 100

-- Define element type
wafer = wafer    or '100'           -- Wafer orientation
angle = angle    or 0               -- Crystal axis angle(radians)
etype = etype    or 'planestrain'   -- Stress analysis type
--post_material = post_material or silicon
--disk_material = disk_material or silicon
post_material = post_material or 'sc_silicon'
disk_material = disk_material or 'sc_silicon'
pelt = mesh:PMLElastic2d_te(post_material)
delt = mesh:PMLELastic2d_te(disk_material)

-- Set default values for any unset parameters
f0    = f0    or 40              -- Stretch function parameter
rpost = rpost or  1e-6           -- Post radius
rdisk = rdisk or 10e-6           -- Disk radius
hdisk = hdisk or  2e-6           -- Height of disk layer
hpost = hpost or  5e-7           -- Height of post layer
rbd   = rbd   or  2e-6           -- Radius of domain
rpml  = rpml  or  4e-6           -- Radius of PML



-- Passing values to matlab
disk_material_n=get_material(disk_material)
disk_material_n=ted_material_normalize(disk_material_n)

-- Mesh consists of three superblocks tied together
mesh:blocks2d( {0, rpost, rbd, rpml},  {-rpml, -rbd, 0}, pelt, order, dense )
mesh:blocks2d( {0, rpost}, {0, hpost, hpost + hdisk},    pelt, order, dense )
mesh:blocks2d( {rpost, rdisk},  {hpost, hpost + hdisk},  delt, order, dense )

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
-- Both disp and thermal bc
function bc_function(x,y)
  if mesheq(x,0)                        then return 'u  ', 0;     end
  if mesheq(x,rdisk)  and meshsg (y,0)  then return 'f  ', 1;     end
  if mesheq(x, rpml)  and meshleq(y,0)  then return 'uuu', 0,0,0; end
  if mesheq(y,-rpml)                    then return 'uuu', 0,0,0; end
  if meshgeq(x,  rbd) and meshleq(y,0)  then return '  u',     0; end
  if meshleq(y, -rbd)                   then return '  u',     0; end 
end
mesh:set_bc(bc_function)
