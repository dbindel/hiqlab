

require 'common.lua'

-- Set default values for any unset paramters

order    = order     or   1              -- Order of shape functions
dense    = dense     or   order          -- Approximate node spacing
f0       = f0        or   20             -- Stretch function
issym    = issym     or    1             -- Whether symmetric anchors
dkRi     = dkRi      or   90e-6
dkWidth  = dkWidth   or   20e-6
bmWidth  = bmWidth   or    6e-6
bmLength = bmLength  or   15e-6
ntDepth  = ntDepth   or   10e-6
ntWidth  = ntWidth   or   10e-6
anLength = anLength  or   25e-6
anWidth  = anWidth   or   50e-6
pmlDepth = pmlDepth  or   10e-6
-- etype    = etype     or   'planestress'  -- Stress analysis type
dense    = dense*dkWidth
dense    = 1e-6
meshtol  = dense/100
bctol    = dense/10
dense_x  = dense
dense_y  = dense
caxis    = { 1, 0, 0, 0, 1, 0}           -- Crystal orientation

if ntWidth < bmWidth then
  ntWidth = bmWidth
end

-- Number of divisions
ndiv_Rim =  ndiv_Rim or 8     -- Through thickness division Rim
ndiv_Rom =  ndiv_Rom or 8     -- Through thickness division Rom
ndiv_Rp  =  ndiv_Rp  or 32    -- Divisions along quarter arc
ndiv_bmW =  3    -- Through bmWidth   division 
ndiv_bmL = 10    -- Along bmLength division
ndiv_ntW =  2    -- Divisions along notch
ndiv_anW = 10    -- Divisions along anchor width/2
ndiv_anL =  8    -- Divisions along anchor depth
ndiv_pml =  8    -- Divisions in pml

-- Nondimensionalize length parameters
--beam_material = beam_material or 'poly_silicon'
beam_material = beam_material or 'sc_silicon'
dim_scales    = ted_nondim(beam_material, dkWidth)

-- To pass values into matlab
beam_material_n = get_material(beam_material)
dims_scales     = ted_material_normalize(beam_material_n)
cL              = dim_scales.L
cT              = dim_scales.T
cM              = dim_scales.M
cTh             = dim_scales.Th

mesh = Mesh:new(2)
belt = mesh:PMLElastic2d_te_planestress(beam_material, caxis)

-- Define mesh function
function open_ring2d(r_in, r_out, slit_l, slit_r, ndiv_r, ndiv_t, elt, order)
  --
  -- mesh open ring with slit at 180deg (and 0 deg)  
  -- slit_l: width of left  slit
  -- slit_r: width of right slit
  -- ndiv_r: divisions along the radius
  -- ndiv_t: divisions along the theta/per quadrant

  local x1, x2, x3, x4
  local ang_l_in, ang_l_out
  local ang_r_in, ang_r_out

  ang_l_in = asin(slit_l/2/r_in)
  ang_l_out= asin(slit_l/2/r_out)
  ang_r_in = asin(slit_r/2/r_in)
  ang_r_out= asin(slit_r/2/r_out)
  x1 = -r_out*cos(ang_l_out)
  x2 = -r_in *cos(ang_l_in )
  x3 =  r_in *cos(ang_r_in )
  x4 =  r_out*cos(ang_r_out)

  mesh:add_curved_block_shape2d(ndiv_r, ndiv_t, elt, order, 
      { x1, slit_l/2, 0, r_out, 
        x2, slit_l/2, 0, r_in }, { 0.0, -1/r_in, 0.0, 1/r_out })
  mesh:add_curved_block_shape2d(ndiv_r, ndiv_t, elt, order, 
      {  0, -r_out, x1, -slit_l/2,
         0, -r_in , x2, -slit_l/2},{ 0.0,-1/r_in, 0.0, 1/r_out })
  mesh:add_curved_block_shape2d(ndiv_r, ndiv_t, elt, order, 
      { x3, slit_r/2, 0, r_in, 
        x4, slit_r/2, 0, r_out}, { 0.0,  1/r_out, 0.0,-1/r_in })
  mesh:add_curved_block_shape2d(ndiv_r, ndiv_t, elt, order, 
      {  0, -r_in , x3, -slit_r/2,
         0, -r_out, x4, -slit_r/2},{ 0.0, 1/r_out, 0.0,-1/r_in })
end

-- Coordinates for Ring
dkRm=  dkRi  + dkWidth - ntDepth
dkRo=  dkRi  + dkWidth
ang_in = asin(ntWidth/2/dkRi)
ang_out= asin(ntWidth/2/dkRm)
x1     = -dkRm*cos(ang_out)
x2     = -dkRi*cos(ang_in)
x3     = -x2
x4     = -x1

-- Coordinates for PML Anchor
xa1 = -anLength - bmLength + x1
xa2 = xa1 + pmlDepth
xa3 = - bmLength + x1
xa4 = -xa3
xa5 = -xa2
xa6 = -xa1

ya1 = - anWidth/2
ya2 = ya1 + pmlDepth
ya3 = -bmWidth/2
ya4 = -ya3
ya5 = -ya2
ya6 = -ya1

-- Mesh inner and outer ring
open_ring2d(dkRi, dkRm, ntWidth, ntWidth*issym, ndiv_Rim, ndiv_Rp, belt, order)
if dkRm~=dkRo then
  open_ring2d(dkRm, dkRo, ntWidth, ntWidth*issym, ndiv_Rom, ndiv_Rp, belt, order)
end

-- Mehs connection Beam
mesh:blocks2dn({ x1-bmLength, x1, x2 }, {1+order*ndiv_bmL, 1+order*ndiv_Rim}, 
               {-bmWidth/2, bmWidth/2}, {1+order*ndiv_bmW}, belt,order)
if issym==1 then
  mesh:blocks2dn({ x3, x4, x4+bmLength }, {1+order*ndiv_Rim, 1+order*ndiv_bmL}, 
                 {-bmWidth/2, bmWidth/2}, {1+order*ndiv_bmW}, belt,order)
end

-- Mesh notched portion if it exists
if ntWidth > bmWidth then
  mesh:blocks2dn({x1, x2}, {ndiv_Rim*order + 1}, 
                 {-ntWidth/2, -bmWidth/2},
                 {1+order*ndiv_ntW},belt, order)
  mesh:blocks2dn({x1, x2}, {ndiv_Rim*order + 1}, 
                 {bmWidth/2, ntWidth/2},
                 {1+order*ndiv_ntW},belt, order)
  if issym==1 then
    mesh:blocks2dn({ x3,  x4}, {ndiv_Rim*order + 1}, 
                   {-ntWidth/2, -bmWidth/2},
                   {1+order*ndiv_ntW},belt, order)
    mesh:blocks2dn({ x3,  x4}, {ndiv_Rim*order + 1}, 
                   {bmWidth/2, ntWidth/2},
                   {1+order*ndiv_ntW},belt, order)
  end
end

-- Mesh anchor
mesh:blocks2dn({xa1, xa2, xa3}, {1+order*ndiv_pml,1+order*ndiv_anL}, 
               {ya1,ya2,ya3,ya4,ya5,ya6},
               {1+order*ndiv_pml,1+order*ndiv_anW,1+order*ndiv_bmW,
                1+order*ndiv_anW,1+order*ndiv_pml},belt,order)
if issym==1 then
  mesh:blocks2dn({xa4, xa5, xa6}, {1+order*ndiv_anL,1+order*ndiv_pml}, 
                 {ya1,ya2,ya3,ya4,ya5,ya6},
                 {1+order*ndiv_pml,1+order*ndiv_anW,1+order*ndiv_bmW,
                  1+order*ndiv_anW,1+order*ndiv_pml},belt,order)
end


-- Mesh tie
mesh:tie()

-- Set stretch function
belt:set_stretch(function(x,y)
  local xs = 0  -- x stretch value
  local ys = 0  -- y stretch value
  if meshleq(x,xa2) then xs = f0*(xa2 - x)/pmlDepth end
  if meshleq(x,xa3) and meshleq(y,ya2) then ys = f0*(ya2 - y)/pmlDepth end
  if meshleq(x,xa3) and meshgeq(y,ya5) then ys = f0*(y - ya5)/pmlDepth end
  -- For symmetric case
  if meshgeq(x,xa5) then xs = f0*(x - xa5)/pmlDepth end
  if meshgeq(x,xa4) and meshleq(y,ya2) then ys = f0*(ya2 - y)/pmlDepth end
  if meshgeq(x,xa4) and meshgeq(y,ya5) then ys = f0*(y - ya5)/pmlDepth end

  return xs,ys
end)

-- Set boundary condition
mesh:set_bc(function(x,y)
  if mesheq( x, xa1) then return 'uuu', 0, 0, 0; end
  if meshleq(x, xa3) and (mesheq( y, ya1) or mesheq(y, ya6))
                     then return 'uuu', 0, 0, 0; end

  if mesheq( x, xa6) then return 'uuu', 0, 0, 0; end
  if meshgeq(x, xa4) and (mesheq( y, ya1) or mesheq(y, ya6))
                     then return 'uuu', 0, 0, 0; end

  -- Fix thermal degrees of freedom in PML
  if f0>0 then
    if meshleq(x, xa2) or meshgeq(x, xa5) then return '  u',  0; end
    if (meshleq(x, xa3) or meshgeq(x, xa4)) and
       (meshleq(y, ya2) or meshgeq(y, ya5))
                       then return '  u',       0; end
  end

  -- Forcing term
  local rad
  rad = sqrt(x^2+y^2)
--  if mesheq(rad, dkRi+dkWidth,bctol) and not
--     (x < 0 and meshbetween(y,-ntWidth/2, ntWidth/2)) then 
--                              return 'ff ', x/rad, y/rad; end
  if mesheq(rad, dkRi+dkWidth,bctol) then return 'ff ', x/rad, y/rad; end
--  if mesheq(rad, dkRi) then return 'ff ', x/rad, y/rad; end
end)
