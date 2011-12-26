

require 'common.lua'

-- Set default values for any unset paramters

order    = order     or   1              -- Order of shape functions
dense    = dense     or   order          -- Approximate node spacing
f0       = f0        or   20             -- Stretch function
dkRo     = dkRo      or    50e-6
bmWidth  = bmWidth   or     5e-6
bmLength = bmLength  or   140e-6
anLength = anLength  or    25e-6
anWidth  = anWidth   or    50e-6
pmlDepth = pmlDepth  or    10e-6
-- etype    = etype     or   '2hd'          -- Stress analysis type
dense    = dense*bmWidth
meshtol  = dense/100
dense_x  = dense
dense_y  = dense
caxis    = { 1, 0, 0, 0, 1, 0}           -- Crystal orientation

-- Mesh input
ndiv_Rom =  ndiv_Rom or  9    -- Through thickness division Rom
ndiv_Rp  =  ndiv_Rp  or 10    -- Divisions along quarter arc
ndiv_bmW =  ndiv_bmW or  3    -- Through bmWidth   division
ndiv_bmL =  ndiv_bmL or 10    -- Along bmLength division
ndiv_anW =  ndiv_anW or  5    -- Divisions along anchor width/2
ndiv_anL =  ndiv_anL or  7    -- Divisions along anchor depth
ndiv_pml =  ndiv_pml or  7    -- Divisions in pml

-- Nondimensionalize length parameters
beam_material = beam_material or 'aln'
dim_scales    = pz_nondim(beam_material, bmWidth)

-- To pass values into matlab
beam_material_n = get_material(beam_material)
dims_scales     = pz_material_normalize(beam_material_n)
cL              = dim_scales.L
cT              = dim_scales.T
cM              = dim_scales.M
cA              = dim_scales.A

-- Construct mesh and material
mesh = Mesh:new(2)
belt = mesh:PMLElastic2hd_pz(beam_material, caxis)

b_len=  sqrt( dkRo^2 - (bmWidth/2)^2 )

-- Top Right
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      { dkRo/2, bmWidth/2, dkRo/2*cos(pi/4), dkRo/2*sin(pi/4),
        b_len , bmWidth/2, dkRo  *cos(pi/4), dkRo  *sin(pi/4)},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      { dkRo/2*cos(pi/4), dkRo/2*sin(pi/4), 0, dkRo/2,
        dkRo  *cos(pi/4), dkRo  *sin(pi/4), 0, dkRo},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rp, ndiv_Rp, belt, order,
      {     0, bmWidth/2,   0,  dkRo/2, 
       dkRo/2, bmWidth/2, dkRo/2*cos(pi/4), dkRo/2*sin(pi/4)})

-- Top Left
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      {     0, dkRo/2, -dkRo/2*cos(pi/4), dkRo/2*sin(pi/4), 
            0, dkRo  , -dkRo  *cos(pi/4), dkRo  *sin(pi/4)},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      { -dkRo/2*cos(pi/4), dkRo/2*sin(pi/4), -dkRo/2, bmWidth/2,
        -dkRo  *cos(pi/4), dkRo  *sin(pi/4), -b_len , bmWidth/2},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rp, ndiv_Rp, belt, order,
      {     0, bmWidth/2, -dkRo/2,  bmWidth/2, 
            0,    dkRo/2, -dkRo/2*cos(pi/4), dkRo/2*sin(pi/4)})

-- Bottom Left
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      {-dkRo/2,-bmWidth/2,-dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4),
       -b_len ,-bmWidth/2,-dkRo  *cos(pi/4),-dkRo  *sin(pi/4)},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      {-dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4), 0,-dkRo/2,
       -dkRo  *cos(pi/4),-dkRo  *sin(pi/4), 0,-dkRo},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rp, ndiv_Rp, belt, order,
      {0, -bmWidth/2,   0, -dkRo/2,
       -dkRo/2,-bmWidth/2, -dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4)})

-- Bottom right
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      {     0,-dkRo/2,  dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4),
            0,-dkRo  ,  dkRo  *cos(pi/4),-dkRo  *sin(pi/4)},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rom, ndiv_Rp, belt, order,
      {  dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4),  dkRo/2,-bmWidth/2,
         dkRo  *cos(pi/4),-dkRo  *sin(pi/4),  b_len ,-bmWidth/2},
      {0, 1/dkRo, 0, 0})
mesh:add_curved_block_shape2d(ndiv_Rp, ndiv_Rp, belt, order,
      {     0,-bmWidth/2,  dkRo/2, -bmWidth/2,
            0,-   dkRo/2,  dkRo/2*cos(pi/4),-dkRo/2*sin(pi/4)})

-- Beam
mesh:blocks2dn({-bmLength/2, -b_len, -dkRo/2, 0, dkRo/2, b_len, bmLength/2}, 
               {ndiv_bmL, ndiv_Rom+1, ndiv_Rp+1,ndiv_Rp+1, ndiv_Rom+1, ndiv_bmL}, 
               {-bmWidth/2, bmWidth/2},
               {ndiv_bmW}, belt,order)

-- PML Anchor
xa1 = - bmLength/2 - anLength
xa2 = xa1 + pmlDepth
xa3 = - bmLength/2
xa4 = -xa3
xa5 = -xa2
xa6 = -xa1
ya1 = - anWidth/2
ya2 = ya1 + pmlDepth
ya3 = -bmWidth/2
ya4 = -ya3
ya5 = -ya2
ya6 = -ya1

mesh:blocks2dn({xa1, xa2, xa3}, {ndiv_pml,ndiv_anL}, {ya1,ya2,ya3,ya4,ya5,ya6},
               {ndiv_pml,ndiv_anW,ndiv_bmW,ndiv_anW,ndiv_pml},belt,order)
mesh:blocks2dn({xa4, xa5, xa6}, {ndiv_pml,ndiv_anL}, {ya1,ya2,ya3,ya4,ya5,ya6},
               {ndiv_pml,ndiv_anW,ndiv_bmW,ndiv_anW,ndiv_pml},belt,order)
mesh:tie()

-- Set stretch function
belt:set_stretch(function(x,y)
  local xs = 0  -- x stretch value
  local ys = 0  -- y stretch value
  if meshleq(x,xa2) then xs = f0*(xa2 - x)/pmlDepth end
  if meshgeq(x,xa5) then xs = f0*(x - xa5)/pmlDepth end
  if meshleq(x,xa3) and meshleq(y,ya2) then ys = f0*(ya2 - y)/pmlDepth end
  if meshleq(x,xa3) and meshgeq(y,ya5) then ys = f0*(y - ya5)/pmlDepth end
  if meshgeq(x,xa4) and meshleq(y,ya2) then ys = f0*(ya2 - y)/pmlDepth end
  if meshgeq(x,xa4) and meshgeq(y,ya5) then ys = f0*(y - ya5)/pmlDepth end

  return xs,ys
end)

-- Set Boundary conditions
mesh:set_bc(function(x,y)

  if mesheq( x, xa1) then return 'uuu', 0, 0, 0; end
  if mesheq( x, xa6) then return 'uuu', 0, 0, 0; end
  if meshleq(x, xa3) and (mesheq( y, ya1) or mesheq(y, ya6))
                     then return 'uuu', 0, 0, 0; end
  if meshgeq(x, xa4) and (mesheq( y, ya1) or mesheq(y, ya6))
                     then return 'uuu', 0, 0, 0; end
  if meshbetween(x, xa1,-b_len) then return '  u',       0; end
  if meshbetween(x, b_len, xa6) then return '  u',       0; end

  if meshbetween(x, 0, b_len) and meshbetween(y, 0, dkRo)
                                then return '  u',       1; end
  if meshbetween(x,-b_len, 0) and meshbetween(y, 0,dkRo) 
                                then return '  u',       0; end
  if meshbetween(x,-b_len, 0) and meshbetween(y, -dkRo, 0) 
                                then return '  u',       0; end
  if meshbetween(x, 0, b_len) and meshbetween(y,-dkRo,0)
                                then return '  u',       0; end
end)

function sensing_vector(x,y)

--  if meshbetween(x, 0, b_len) and meshbetween(y, 0, dkRo)
--                                then return '  u',       0; end
--  if meshbetween(x,-b_len, 0) and meshbetween(y, 0,dkRo) 
--                                then return '  u',       0; end
  if meshbetween(x,-b_len, 0) and meshbetween(y, -dkRo, 0) 
                                then return '  u',       1; end
--  if meshbetween(x, 0, b_len) and meshbetween(y,-dkRo,0)
--                                then return '  u',       0; end
end
