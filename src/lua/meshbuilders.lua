
--
-- Set x -> tfunc(x) for all nodes between nstart and nend
--
function Mesh:transform_nodes(nstart,nend,tfunc)
  for j = nstart,nend do
    self:set_x(j,{tfunc(self:x(j))})
  end
end


--
-- Rewrite node i -> node tfunc(i) in ix array for elts nstart to nend
-- Only affects the first nen elements of the IX array
--
function Mesh:transform_conn(nstart,nend,nen,tfunc)
  for j = nstart,nend do
    local ix_old = {}
    local ix_new = {}
    for i = 1,nen do
      ix_old[i] = self:ix(i-1,j)
    end
    ix_new = tfunc(ix_old)
    for i = 1,nen do
      self:set_ix(i-1,j,ix_new[i])
    end
  end
end


--
-- Add a block on [-1,1]^ndm and transform
--
function Mesh:add_block_transform(...)
  local f = arg[arg.n]
  local nstart = self:numnp()
  arg.n = arg.n - 1
  if arg.n == 4 then
    self:add_block(-1,-1,    1,1,   unpack(arg))
  else
    self:add_block(-1,-1,-1, 1,1,1, unpack(arg))
  end
  local nend = self:numnp()
  self:transform_nodes(nstart,nend-1,f)
end


--
-- Add a block on [-1,1]^ndm and transform by coordinate interpolation
--
function Mesh:add_block_shape(...)
  local f = arg[arg.n]
  arg[arg.n] = function(...)
    interp_coord(f,arg)
    return unpack(arg)
  end
  self:add_block_transform(unpack(arg))
end


local function arc_interpolate(x1, x2, curv, num_node, x)

  local ndiv = num_node-1
  local x_ind= floor((x+1)/2*ndiv)
  local d    = sqrt( (x1[1]-x2[1])^2 + (x1[2]-x2[2])^2 )
  local beta = 2*asin(0.5*d*curv)/ndiv

  local sb, cb = 0, 0
  for i = 0,ndiv-1 do
    cb = cb + cos(beta*i)
    sb = sb + sin(beta*i)
  end
  local det = cb^2 + sb^2
  local vec = {(cb*(x2[1]-x1[1])+sb*(x2[2]-x1[2]))/det,
              (-sb*(x2[1]-x1[1])+cb*(x2[2]-x1[2]))/det}

  sb, cb = 0, 0
  for i = 0,x_ind-1 do
    cb = cb + cos(beta*i)
    sb = sb + sin(beta*i)
  end

  return x1[1] + cb*vec[1]-sb*vec[2],
         x1[2] + sb*vec[1]+cb*vec[2]
end


--
-- Add a block on [-1,1]^ndm with given curvature on each edge.
--
function Mesh:add_curved_block_shape2d(nx, ny, material, order1, xn, curv)
 self:add_block_transform(nx, ny, material, order1,
   function(x,y)
     local a1 = {arc_interpolate({xn[1],xn[2]}, {xn[3],xn[4]},-curv[1], ny, y)}
     local a2 = {arc_interpolate({xn[5],xn[6]}, {xn[7],xn[8]}, curv[2], ny, y)}
     local cv = curv[3] + (-curv[4]-curv[3])*(y+1)/2
     return arc_interpolate(a1, a2, cv, nx, x)
   end)
end


-- "Superblock" generator with pre-set node spacings

function Mesh:blocks1d(xlist, material, order1, dense1)
  order1 = order1 or order
  dense1 = dense1 or dense
  for i = 1,table.getn(xlist)-1 do
    self:add_block( xlist[i  ], xlist[i+1],
                    1 + order1*math.ceil((xlist[i+1]-xlist[i]) / dense1),
                    material, order1 )
  end
end

function Mesh:blocks2d(xlist, ylist, material, order1, dense1, dense2)
  order1 = order1 or order
  dense1 = dense1 or dense
  dense2 = dense2 or dense1
  for i = 1,table.getn(xlist)-1 do
    for j = 1,table.getn(ylist)-1 do
      self:add_block( xlist[i  ], ylist[j  ],
                      xlist[i+1], ylist[j+1],
                      1 + order1*math.ceil((xlist[i+1]-xlist[i]) / dense1),
                      1 + order1*math.ceil((ylist[j+1]-ylist[j]) / dense2),
                      material, order1 )
    end
  end
end

function Mesh:blocks3d(xlist, ylist, zlist, material, order1, dense1,
                                                      dense2, dense3)
  order1 = order1 or order
  dense1 = dense1 or dense
  dense2 = dense2 or dense1
  dense3 = dense3 or dense2
  for i = 1,table.getn(xlist)-1 do
    for j = 1,table.getn(ylist)-1 do
      for k = 1,table.getn(zlist)-1 do
          self:add_block( xlist[i  ], ylist[j  ], zlist[k  ],
                          xlist[i+1], ylist[j+1], zlist[k+1],
                          1 + order1*math.ceil((xlist[i+1]-xlist[i]) / dense1),
                          1 + order1*math.ceil((ylist[j+1]-ylist[j]) / dense2),
                          1 + order1*math.ceil((zlist[k+1]-zlist[k]) / dense3),
                          material, order1 )
      end
    end
  end
end

-- "Superblock" generator with number of nodes specified

function Mesh:blocks1dn(xlist, mlist, material, order1)
  order1 = order1 or order
  for i = 1,table.getn(xlist)-1 do
    self:add_block( xlist[i  ], xlist[i+1],
                    mlist[i],
                    material, order1 )
  end
end


function Mesh:blocks2dn(xlist, mlist, ylist, nlist, material, order1)
  order1 = order1 or order
  for i = 1,table.getn(xlist)-1 do
    for j = 1,table.getn(ylist)-1 do
      self:add_block( xlist[i  ], ylist[j  ],
                      xlist[i+1], ylist[j+1],
                      mlist[i], nlist[j],
                      material, order1 )
    end
  end
end

function Mesh:blocks3dn(xlist, mlist, ylist, nlist, zlist, plist, 
                        material, order1)
  order1 = order1 or order
  for i = 1,table.getn(xlist)-1 do
    for j = 1,table.getn(ylist)-1 do
      for k = 1,table.getn(zlist)-1 do
         self:add_block( xlist[i  ], ylist[j  ], zlist[k  ],
                         xlist[i+1], ylist[j+1], zlist[k+1],
                         mlist[i], nlist[j], plist[k],
                         material, order1 )
      end
    end
  end
end


