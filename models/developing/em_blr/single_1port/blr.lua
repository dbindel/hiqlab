-- Construct middle strip in local coordinates
function middle_strip()

    -- Location of dielectric strips
    local arrd = {-Dx-Dt,-Dx,Dx,Dx+Dt}

    -- Nodal control points in direction x
    local arrx = {-Bw/2,Bw/2}

    -- Nodal control points in direction y
    local arry = {-Bh/2,-Bh/2+Cw,-Cw/2,Cw/2,Bh/2-Cw,Bh/2,
                  -Dx-Dt,-Dx,Dx,Dx+Dt}
    table.sort(arry)
    arry = table.pack_numbers(arry,meshtol or 0)

    -- Find indices corresponding to dielectric position
    local dind1,dind2
    for i=1,table.getn(arry) do
        if mesheq(arry[i],arrd[1]) then
           dind1 = i
        end
        if mesheq(arry[i],arrd[3]) then
           dind2 = i
        end
    end

    -- Mesh the middle strip
    local arry1 = {}
    local arry2 = {}
    local arry3 = {}
    for i=1,dind1 do
        arry1[i] = arry[i]
    end
    for i=dind1+1,dind2 do
        arry2[i-dind1] = arry[i]
    end
    for i=dind2+1,table.getn(arry) do
        arry3[i-dind2] = arry[i]
    end

    -- Construct mesh
    mesh:blocks2d(arrx,            arry1,etype_c, order,densex,densey)
    mesh:blocks2d(arrx,{arrd[1],arrd[2]},etype_em,order,densex,densey)
    mesh:blocks2d(arrx,            arry2,etype_c, order,densex,densey)
    mesh:blocks2d(arrx,{arrd[3],arrd[4]},etype_em,order,densex,densey)
    mesh:blocks2d(arrx,            arry3,etype_c, order,densex,densey)

    -- Fill gap with elastic material
    if fill_gap==1 then
        mesh:blocks2d(arrx,{arrd[1],arrd[2]},etype_g ,order,densex,densey)
        mesh:blocks2d(arrx,{arrd[3],arrd[4]},etype_g ,order,densex,densey)
    end
end

-- Construct supporting beam in local coordinates
function beam_connection(btype,xn,yn)

    local alpha_b= find_angle(xn[3]-xn[1],yn[3]-yn[1])
    local xbc= (xn[1]+xn[4])/2
    local ybc= (yn[1]+yn[4])/2
    local bl = sqrt( (xn[3]-xn[1])^2 + (yn[3]-yn[1])^2 )
    local bw = sqrt( (xn[2]-xn[1])^2 + (yn[2]-yn[1])^2 )

    local x1,x2,x3,x4,x5,x6
    local y1,y2,y3,y4
    x1 = -bl/2
    x2 = -bl/6-bw/2
    x3 = -bl/6+bw/2
    x4 = -x3
    x5 = -x2
    x6 = -x1
    y1 = -bw/2
    y2 =  bw/2
    y3 =  bl/3-bw/2
    y4 =  bl/3+bw/2

    local nstart = mesh:numnp()
    if btype=='straight' then
       mesh:blocks2d({x1,            x6},{y1,y2      },etype_c,order,dense)
    elseif btype=='serpentine' then
       mesh:blocks2d({x1,x2,x3         },{y1,y2      },etype_c,order,dense)
       mesh:blocks2d({   x2,x3,x4,x5   },{      y3,y4},etype_c,order,dense)
       mesh:blocks2d({         x4,x5,x6},{y1,y2      },etype_c,order,dense)
       mesh:blocks2d({   x2,x3         },{   y2,y3   },etype_c,order,dense)
       mesh:blocks2d({         x4,x5   },{   y2,y3   },etype_c,order,dense)
    end
    local nend = mesh:numnp()

    local function transform_b(x,y)
        local x1, y1
        x1 = xbc+ x*cos(alpha_b) - y*sin(alpha_b)
        y1 = ybc+ x*sin(alpha_b) + y*cos(alpha_b)
        return x1,y1
    end

    -- Transform the coordinates
    mesh:transform_nodes(nstart,nend-1,transform_b)

end

-- Construct legs in local coordinates
function construct_legs()

    -- Set arrays
    local x1,x2,x3,x4
    local y1,y2,y3,y4,y5,y6
    x1 = -Bw/2-Cl
    x2 = -Bw/2
    x3 =  Bw/2
    x4 =  Bw/2+Cl
    y1 = -Bh/2
    y2 = -Bh/2+Cw
    y3 =  -Cw/2
    y4 =   Cw/2
    y5 =  Bh/2-Cw
    y6 =  Bh/2

    local xlegs = {x1,x1,x2,x2,
                   x3,x3,x4,x4}
    local ylegs = {y3,y4,y3,y4,
                   y3,y4,y3,y4,
                   y5,y6,y5,y6,
                   y5,y6,y5,y6,
                   y1,y2,y1,y2,
                   y1,y2,y1,y2}

    -- Construct legs
    for i = 1,6 do
        local j = mod(i+1,2)*4
        local k = 4*(i-1)
        beam_connection(Csh,{xlegs[j+1],xlegs[j+2],xlegs[j+3],xlegs[j+4]},
                            {ylegs[k+1],ylegs[k+2],ylegs[k+3],ylegs[k+4]})
    end

end

-- Construct anchors in local coordinates
function construct_anchors()

    -- Return values
    local bc_t = {} 
    local st_t = {}

    -- Set arrays
    local x1,x2,x3,x4
    local y1,y2,y3,y4,y5,y6
    x1 = -Bw/2-Cl
    x2 = -Bw/2
    x3 =  Bw/2
    x4 =  Bw/2+Cl
    y1 = -Bh/2
    y2 = -Bh/2+Cw
    y3 =  -Cw/2
    y4 =   Cw/2
    y5 =  Bh/2-Cw
    y6 =  Bh/2

    local xanchor={x1,
                   x4}
    local yanchor={y4,y3,
                   y3,y4,
                   y2,y1,
                   y1,y2,
                   y6,y5,
                   y5,y6}

    for i=1,6 do
        local j = mod(i+1,2)+1
        local k = 2*(i-1)
        local x1,y1 = xanchor[j],yanchor[k+1]
        local x2,y2 = xanchor[j],yanchor[k+2]
        bc_t[i],st_t[i] = 
               pml_blocks2d({x1,x2},{y1,y2},4,Pw,Ph,Pd,f0,etype_c,order,dense,dense)
    end

    return bc_t,st_t
end

-- Construct global plate shape functions in local coordinates
function construct_1port_plate_functions()

    local gs_t = {} 
    local xg1,xg2,yg1,yg2,yg3,yg4
    xg1 = -Bw/2    
    xg2 =  Bw/2   
    yg1 = -Dx - Dt 
    yg2 = -Dx     
    yg3 =  Dx      
    yg4 =  Dx + Dt 

    -- Drive electrode
    gs_t[1] = function(x,y)
                  if mesheq(y, yg1) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
                  if mesheq(y, yg4) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
              end
    -- Sense electrode
    gs_t[2] = function(x,y)
                  if mesheq(y, yg2) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
                  if mesheq(y, yg3) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
              end

    return gs_t
end

function construct_2port_plate_functions()

    local gs_t = {} 
    local xg1,xg2,yg1,yg2,yg3,yg4
    xg1 = -Bw/2    
    xg2 =  Bw/2   
    yg1 = -Dx - Dt 
    yg2 = -Dx     
    yg3 =  Dx      
    yg4 =  Dx + Dt 

    -- Drive electrode
    gs_t[1] = function(x,y)
                  if mesheq(y, yg1) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
              end

    -- Sense electrode
    gs_t[2] = function(x,y)
                  if mesheq(y, yg4) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
              end

    -- Base electrode
    gs_t[3] = function(x,y)
                  if mesheq(y, yg2) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
                  if mesheq(y, yg3) and meshbetween(x, xg1, xg2)
                                    then return    0,0,1,0; end
              end

    return gs_t
end


-- Construct bulk lateral resonator in global coordinates     
function construct_blr(port,xc,yc,alpha)

    -- Return values
    local bc_func = {}
    local st_func = {}
    local gs_func = {}

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

    -- Start construction of mesh
    local nstart = mesh:numnp()

    -- Construct middle strip in local coordinates
    middle_strip()

    -- Construct legs in local coordinates
    construct_legs()

    -- Construct anchors functions in local coordinates
    local bc_t,st_t = construct_anchors()

    local nend = mesh:numnp()

    -- Transform the mesh to global coordinates
    mesh:transform_nodes(nstart,nend-1,transform)

    -- Construct global shape functions in local coordinates
    local gs_t
    if port == 1 then
        gs_t = construct_1port_plate_functions()
    elseif port == 2 then
        gs_t = construct_2port_plate_functions()
    end

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

    for i=1,table.getn(gs_t) do
        local gs = gs_t[i]
        gs_func[i] = function(x,y)
                         local x3,y3 = inv_transform(x,y)
                         return gs(x3,y3)
                     end
    end

    return bc_func,st_func,gs_func
end

