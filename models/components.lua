-- Functions for inverting and reflecting nodes
function move_node(x1,y1,rfpt)
  -- rfpt: new origin 
    rfpt = rfpt or { 0, 0}
    return x1+rfpt[1],y1+rfpt[2]
end

function move_nodes(nodex,rfpt)
    for i = 1,table.getn(x),2 do
        newnodex[i], newnodex[i+1] = 
            move_node(nodex[i],nodex[i+1],rfpt)
    end
    return newnodex
end

function rotate_node(x1,y1,theta,rfpt)
  -- rfpt: rotation orgin
    rfpt = rfpt or { 0, 0}
    local xnew = rfpt[1] + cos(theta)*x1 -sin(theta)*y1
    local ynew = rfpt[2] + sin(theta)*x1 +cos(theta)*y1
    
    return xnew,ynew
end

function rotate_nodes(nodex,theta,rfpt)
    local newnodex = {}
    for i = 1,table.getn(nodex),2 do
        newnodex[i], newnodex[i+1] = 
            rotate_node(nodex[i],nodex[i+1],theta,rfpt)
    end
    return newnodex
end

function invert_node(x1,y1,rfpt)
  -- rfpt: reference point for the inversion
    rfpt = rfpt or { 0, 0}
    local xnew = rfpt[0] - (x1 - rfpt[0])
    local ynew = rfpt[1] - (y1 - rfpt[1])

    return xnew,ynew
end

function invert_nodes(nodex, rfpt)
  -- rfpt: reference point for the inversion
    local newnodex = {}
    for i = 1,table.getn(x),2 do
        newnodex[i], newnodex[i+1] = 
            invert_node(nodex[i],nodex[i+1],rfpt)
    end
    return newnodex
end

function reflect_node(x1,y1,nvec,rfpt)
  -- nvec: table of normal vector of plane
  -- rfpt: reference point for the reflection plane
  local rfpt = rfpt or { 0, 0 }
  local norm = sqrt(nvec[1]^2 + nvec[2]^2)
        nvec = {nvec[1]/norm, nvec[2]/norm}


  local t    = nvec[1]*(x1 - rfpt[1]) +
               nvec[2]*(y1 - rfpt[2])

  local ref  = { rfpt[1] + nvec[1]*t,
                 rfpt[2] + nvec[2]*t}
  local xnew = x1 - 2*t*nvec[1]
  local ynew = y1 - 2*t*nvec[2]

  return xnew,ynew
end

function reflect_nodes(nodex, nvec, rfpt)
  -- nvec: table of normal vector of plane
  -- rfpt: reference point for the reflection plane
    local newnodex = {}
    for i = 1,table.getn(nodex),2 do
        newnodex[i], newnodex[i+1] = 
            reflect_node(nodex[i],nodex[i+1],nvec,rfpt)
    end
    return newnodex
end

function reflect_quad(nodex,nvec,rfpt)
  local newx = reflect_nodes(nodex,nvec,rfpt)
--  newx = {newx[3],newx[4],
--          newx[1],newx[2],
--          newx[7],newx[8],
--          newx[5],newx[6]}
  newx = {newx[5],newx[6],
          newx[7],newx[8],
          newx[1],newx[2],
          newx[3],newx[4]}

  return newx
end

-- Functions for generating PML blocks(bc and stretch functions)
function find_angle(x,y)
    local r, alpha
    r     = math.sqrt(x^2 + y^2)
    alpha = math.acos(x/r)
    if y < 0 then
      alpha = -alpha
    end
    return alpha
end

-- Functions for generating PML blocks(bc and stretch functions)
function pml_blocks2d(xnodes,ynodes,ndf,bw,bh,dpml,f0,etype,order,densex,densey)

    local xarr = {}
    local yarr = {}
    local nx   = table.getn(xnodes)
    local ny   = table.getn(ynodes)
    local a1,a2,a3
    local meshtol = meshtol or bw*1e-6
    local bcarray

    if ndf == 2 then
        bcarray = {'uu',0,0}
    elseif ndf == 3 then
        bcarray = {'uuu',0,0,0}
    elseif ndf == 4 then
        bcarray = {'uuuu',0,0,0,0}
    end

    local xc,yc,alpha
    local nx,ny
    nx = table.getn(xnodes)
    ny = table.getn(ynodes)
    xc = (xnodes[1] + xnodes[nx])/2
    yc = (ynodes[1] + ynodes[ny])/2
    alpha = find_angle(xnodes[nx]-xnodes[1],ynodes[ny]-ynodes[1])

    xarr[1   ] = -bw/2
    xarr[2   ] = -bw/2 + dpml
    xarr[nx+3] =  bw/2 - dpml
    xarr[nx+4] =  bw/2
    xarr[3   ] = - sqrt( (xnodes[1]-xc)^2 + (ynodes[1]-yc)^2 )
    for i=2,nx do
        xarr[i+2] = xarr[3] + sqrt( (xnodes[i]-xnodes[1])^2 + (ynodes[i]-ynodes[1])^2 )
    end
    yarr = {-bh,-bh+dpml,0.0}

    -- Define forward and inverse functions
    local function transform(x,y)
        local x1, y1
        x1 = xc + x*cos(alpha) - y*sin(alpha)
        y1 = yc + x*sin(alpha) + y*cos(alpha)
        return x1,y1
    end
    local function inv_transform(x,y)
        local x1, y1
        x1 =  cos(-alpha) * (x-xc) - sin(-alpha) * (y-yc)
        y1 =  sin(-alpha) * (x-xc) + cos(-alpha) * (y-yc)
        return x1,y1
    end

    -- Add a block
    local nstart = mesh:numnp()
    mesh:blocks2d(xarr,yarr,etype,order,densex,densey)
    local nend   = mesh:numnp()
    -- Transform the coordinates
    mesh:transform_nodes(nstart,nend-1,transform)

    -- Define boundary condition function
    local function bc_func(x,y)

        local x1,y1 = inv_transform(x,y)

        if meshbetween(y1,yarr[1],yarr[3],meshtol) and 
          (mesheq(x1,xarr[1],meshtol) or mesheq(x1,xarr[nx+4],meshtol))
            then return unpack(bcarray); end
        if meshbetween(x1,xarr[1],xarr[nx+4],meshtol) and mesheq(y1,yarr[1],meshtol) 
            then return unpack(bcarray); end

    end

    local function s_func(x,y)
        local sx = 0
        local sy = 0

        local x1,y1 = inv_transform(x,y)

        if meshbetween(x1,xarr[1],xarr[2],meshtol) and 
           meshbetween(y1,yarr[1],yarr[3],meshtol) then
           sx = f0*abs(x1-xarr[2])/dpml; end
        if meshbetween(x1,xarr[nx+3],xarr[nx+4],meshtol) and 
           meshbetween(y1,yarr[1],yarr[3],meshtol) then
           sx = f0*abs(x1-xarr[nx+3])/dpml; end
        if meshbetween(x1,xarr[1],xarr[nx+4],meshtol) and 
           meshbetween(y1,yarr[1],yarr[2],meshtol) then
           sy = f0*abs(y1-yarr[2])/dpml; end

        local sx1,sx2 = rotate_node(sx, 0,alpha)
        local sy1,sy2 = rotate_node( 0,sy,alpha)
        sx = abs(sx1) + abs(sy1)
        sy = abs(sx2) + abs(sy2)

        return sx,sy
    end

    return bc_func, s_func
end

function len2div(len,order1,dense1)
    order1 = order1 or order
    dense1 = dense1 or dense
    return max(order1*ceil((len)/dense1),1)
end

function table.pack_numbers(x,tol)
  -- packs a numerical table x depending on tol
  tol = tol or 0
  local newt = {}
  local nx  = table.getn(x)

  if nx==1 then return
  else newt[1] = x[1]
  end

  local ind = 1
  for i = 2,nx do
      if(abs(newt[ind]-x[i])>tol) then
        ind = ind + 1
        newt[ind] = x[i]
      end
  end

  return newt
end
