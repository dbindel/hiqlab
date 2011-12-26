-- Generate Free Free beam in global coordinates
function generate_frfr_beam2d(xc,yc,alpha)

  local x1,x2,x3,x4,x5,x6
  local y1,y2,y3,y4
  x1 = - Bl/2
  x2 = - Sx - Sw
  x3 = - Sx 
  x4 =   Sx 
  x5 =   Sx + Sw
  x6 =   Bl/2
  y1 =  -Sl - Bw/2
  y2 =  -Bw/2
  y3 =   Bw/2
  y4 =   Sl + Bw/2
  bc_func = {}
  st_func = {}

  -- Define transformation
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

  local nstart = mesh:numnp()
  mesh:blocks2d({x1,x2,x3,x4,x5,x6},{   y2,y3   },etype,order,dense)
  mesh:blocks2d({   x2,x3         },{y1,y2      },etype,order,dense)
  mesh:blocks2d({   x2,x3         },{      y3,y4},etype,order,dense)
  mesh:blocks2d({         x4,x5   },{y1,y2      },etype,order,dense)
  mesh:blocks2d({         x4,x5   },{      y3,y4},etype,order,dense)
  local nend   = mesh:numnp()

  -- Transform the coordinates
  mesh:transform_nodes(nstart,nend-1,transform)

  -- Construct anchors
  local xbc = {x2,x3,y1,y1,
               x3,x2,y4,y4,
               x4,x5,y1,y1,
               x5,x4,y4,y4}

  local bc_t = {}
  local st_t = {}
  bc_t[1],st_t[1] = pml_blocks2d({x2,x3},{y1,y1},3,Pw,Ph,Pd,f0,etype,order,dense)
  bc_t[2],st_t[2] = pml_blocks2d({x3,x2},{y4,y4},3,Pw,Ph,Pd,f0,etype,order,dense)
  bc_t[3],st_t[3] = pml_blocks2d({x4,x5},{y1,y1},3,Pw,Ph,Pd,f0,etype,order,dense)
  bc_t[4],st_t[4] = pml_blocks2d({x5,x4},{y4,y4},3,Pw,Ph,Pd,f0,etype,order,dense)

  -- Transform to global coordinate system
  for i=1,table.getn(bc_t) do
      local bc = bc_t[i]
      local st = st_t[i]
      bc_func[i] = function(x,y)
                       local x3,y3 = inv_transform(x,y)
                       return bc(x3,y3)
                   end
      st_func[i] = function(x,y)
                       local x3,y3 = inv_transform(x,y)
                       local sx,sy = st(x3,y3)

                       local sx1,sx2 = rotate_node(sx, 0,alpha)
                       local sy1,sy2 = rotate_node( 0,sy,alpha)
                       sx = abs(sx1) + abs(sy1)
                       sy = abs(sx2) + abs(sy2)

                       return sx,sy
                   end
  end

  return bc_func,st_func
end
