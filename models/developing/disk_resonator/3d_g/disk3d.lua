-- Include function definition file
require 'common.lua'
require 'disk3d_func.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(3);

-- Define mesh tolerance
order   = order or 1               -- Order of shape functions
dense   = dense or order*0.5e-6    -- Approximate element size
meshtol = dense/100                -- Default mesh tolerance
nd      = nd    or 5               -- Number of elements along pi/4

-- Define element type
post_material = post_material or 'silicon2'
disk_material = disk_material or 'silicon2'
etype= mesh:PMLElastic3d(disk_material)

-- Define geometric parameters
f0   = f0 or 40           -- Stretch function parameter
rpost= rpost or  1e-6     -- Post radius
rdisk= rdisk or 10e-6     -- Disk radius
hdisk= hdisk or 2e-6      -- Height of disk layer
hpost= hpost or 5e-7      -- Height of post layer
rbd  = rbd   or 2e-6      -- Radius of domain
rpml = rpml  or 4e-6      -- Radius of pml
xo   = xo    or 0         -- Shift of post
yo   = yo    or 0         -- Shift of post

rt1  = 4e-6;
rt2  = rt1+dense;
rt3  = rt1+dense*2;
rb1  = 1.5e-6;

-- Define mesh using block command
for i = 1, 4 do

ang = pi/2*(i-1)

-- Disk portion
mesh_inner_disk(    xo, yo, rpost,      {hpost, hpost+hdisk}, nd  , ang)
mesh_post_interface(xo, yo, rpost, rt1, {hpost, hpost+hdisk}, nd  , ang)
mesh_transition(rt1, rt2, rt3,          {hpost, hpost+hdisk}, nd  , ang)
mesh_steady(rt3, rdisk,                 {hpost, hpost+hdisk}, nd*3, ang)

-- Stem portion
mesh_inner_disk(xo, yo, rpost, {0, hpost}, nd, ang)

-- Base
mesh_inner_disk(    xo, yo, rpost,      {-rbd, 0}, nd, ang)
mesh_post_interface(xo, yo, rpost, rb1, {-rbd, 0}, nd, ang)
mesh_steady_end(rb1, rbd,               {-rbd, 0}, nd, ang)

-- Anchor
mesh_inner_disk(    xo, yo, rpost,      {-rpml, -rbd}, nd, ang)
mesh_post_interface(xo, yo, rpost, rb1, {-rpml, -rbd}, nd, ang)
mesh_steady_end(rb1, rbd,               {-rpml, -rbd}, nd, ang)

end

mesh_pmls(rbd,rpml)

-- Tie mesh together
mesh:tie()

-- Define stretching function
function stretch_function(x,y,z)

  local xs = 0   -- x stretch value
  local ys = 0   -- y stretch value
  local zs = 0   -- y stretch value

  if abs(x) > rbd and meshleq(z,0) then xs = f0*(x-rbd)/(rpml-rbd) end
  if abs(y) > rbd and meshleq(z,0) then ys = f0*(y-rbd)/(rpml-rbd) end
  if meshleq(z, -rbd)              then zs =-f0*(z+rbd)/(rpml-rbd) end

  return xs, ys, zs
end
etype:set_stretch(stretch_function)

-- Define boundary conditions
function bc_function(x,y,z)
  if mesheq(abs(x), rpml) and meshleq(z,0)  then return 'uuu', 0, 0, 0; end
  if mesheq(abs(y), rpml) and meshleq(z,0)  then return 'uuu', 0, 0, 0; end
  if mesheq(-z    ,-rpml)                   then return 'uuu', 0, 0, 0; end
--  if mesheq(z    , 0)                   then return 'uuu', 0, 0, 0; end
end
mesh:set_bc(bc_function)

--point_load(0,1, {uy = -5})
--clamp_boundary(function(x,y) return mesheq(x,0) end, 'ux')
--clamp_boundary(function(x,y) return mesheq(y,0) end, 'uy')
