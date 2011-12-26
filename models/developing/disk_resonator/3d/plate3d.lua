-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(3);

-- Define mesh tolerance
order   = order or 1               -- Order of shape functions
dense   = dense or order*0.5e-6    -- Approximate element size
meshtol = dense/100                -- Default mesh tolerance
nd      = nd    or 5               -- Number of elements along pi/4
densez  = densez or dense          -- Approximate elemnent size in PML regions

-- Define element type
post_material = post_material or 'silicon2'
disk_material = disk_material or 'silicon2'
etype= mesh:PMLElastic3d(disk_material)

-- Define geometric parameters
f0   = f0 or 40           -- Stretch function parameter
rpost= rpost or  1e-0     -- Post radius
rdisk= rdisk or 10e-0     -- Disk radius
hdisk= hdisk or 2e-0      -- Height of disk layer
hpost= hpost or 5e-0      -- Height of post layer
rbd  = rbd   or 2e-0      -- Radius of domain
rpml = rpml  or 4e-0      -- Radius of pml
xo   = xo    or 0         -- Shift of post
yo   = yo    or 0         -- Shift of post

-- Define mesh using block command
    local ndivplt = ceil((rdisk-rpost)/dense)*order+1
    local ndivpst = ceil( rpost       /dense)*order+1
    local ndivptk = ceil( hdisk       /dense)*order+1
    local ndivpht = ceil( hpost       /dense)*order+1
    local ndivbet = ceil((rbd-rpost)  /dense)*order+1
    local ndivpml = ceil((rpml-rbd)/densez)*order+1
    local ndivbet = ceil( rbd         /dense)*order+1

--mesh:blocks3dn( {   -rdisk, -rpost, 0, rpost, rdisk}, { ndivplt, ndivpst, ndivpst, ndivplt }, 
--                {   -rdisk, -rpost, 0, rpost, rdisk}, { ndivplt, ndivpst, ndivpst, ndivplt }, 
--                {    hpost,  hpost+hdisk           }, { ndivptk}, etype, order)
mesh:blocks3dn( {           -rpost, 0, rpost       }, {          ndivpst, ndivpst          }, 
                {           -rpost, 0, rpost       }, {          ndivpst, ndivpst          }, 
                {        0,  hpost                 }, { ndivpht}, etype, order)
mesh:blocks3dn( {   -rpml, -rbd, -rpost, 0, rpost, rbd, rpml     }, {ndivpml, ndivbet, ndivpst, ndivpst, ndivbet, ndivpml}, 
                {   -rpml, -rbd, -rpost, 0, rpost, rbd, rpml     }, {ndivpml, ndivbet, ndivpst, ndivpst, ndivbet, ndivpml},
                {   -rpml, -rbd, 0}, {ndivpml, ndivbet}, etype, order)

-- Tie mesh together
mesh:tie()

-- Define stretching function
function stretch_function(x,y,z)

  local xs = 0   -- x stretch value
  local ys = 0   -- y stretch value
  local zs = 0   -- y stretch value

  if abs(x) > rbd and meshleq(z,0) then xs = f0*(abs(x)-rbd)/(rpml-rbd) end
  if abs(y) > rbd and meshleq(z,0) then ys = f0*(abs(y)-rbd)/(rpml-rbd) end
  if meshleq(z, -rbd)              then zs = f0*(abs(z)-rbd)/(rpml-rbd) end

  return xs, ys, zs
end
--etype:set_stretch(stretch_function)

-- Define boundary conditions
function bc_function(x,y,z)
  if mesheq(abs(x), rpml) and meshleq(z,0)  then return 'uuu', 0, 0, 0; end
  if mesheq(abs(y), rpml) and meshleq(z,0)  then return 'uuu', 0, 0, 0; end
  if mesheq( z    ,-rpml)                   then return 'uuu', 0, 0, 0; end
  if mesheq( z    ,hdisk+hpost)             then return 'f  ', 1;       end
end
mesh:set_bc(bc_function)

--point_load(0,1, {uy = -5})
--clamp_boundary(function(x,y) return mesheq(x,0) end, 'ux')
--clamp_boundary(function(x,y) return mesheq(y,0) end, 'uy')
