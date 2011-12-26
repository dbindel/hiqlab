-- Include function definition file
require 'common.lua'

function mesh_steady_end1(ri,ro,zt,nd,angi,ipart)

    local angi   = angi or 0
    local ipart  = ipart or 1
    local numnr  = ceil((ro-ri)/dense)*order+1
    local numns  = order+1
    local numnz  = ceil((zt[2]-zt[1])/dense)*order+1
    local rix = {}
    local riy = {}
    local rox = {}
    local roy = {}
    local xc = {}

    if ipart==1 then
        xc[1] = ro*cos(angi)
        xc[2] = ro*sin(angi)
        xc[3] = xc[1] + ro*cos(angi+pi/2) 
        xc[4] = xc[2] + ro*sin(angi+pi/2) 
    elseif ipart==2 then
        xc[3] = ro*cos(angi+pi/2)
        xc[4] = ro*sin(angi+pi/2)
        xc[1] = xc[3] + ro*cos(angi) 
        xc[2] = xc[4] + ro*sin(angi)
        angi  = angi + pi/4
    end

    -- Nodes in inner side
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        rix[i] = ri*cos(ang)
        riy[i] = ri*sin(ang)
    end

    -- Nodes in outer side
    for i = 1, nd+1 do
        rox[i] = (nd+1 - i)/nd * xc[1] + (i-1)/nd * xc[3]
        roy[i] = (nd+1 - i)/nd * xc[2] + (i-1)/nd * xc[4]
    end

    -- Construct mesh
    for i = 1, nd do
        local xc = {rix[i],   riy[i],   zt[1],
                    rix[i],   riy[i],   zt[2],
                    rix[i+1], riy[i+1], zt[1],
                    rix[i+1], riy[i+1], zt[2],
                    rox[i],   roy[i],   zt[1],
                    rox[i],   roy[i],   zt[2],
                    rox[i+1], roy[i+1], zt[1],
                    rox[i+1], roy[i+1], zt[2]}
        mesh:add_block_shape(numnr, numns, numnz, etype, order, xc)
    end

end

function mesh_steady_end(ri,ro,zt,nd,angi)
   mesh_steady_end1(ri,ro,zt,nd,angi,1)
   mesh_steady_end1(ri,ro,zt,nd,angi,2)
end

function mesh_steady1(ri,ro,zt,nd,angi)

    local angi   = angi or 0
    local numnr  = ceil((ro-ri)/dense)*order+1
    local numns  = order+1
    local numnz  = ceil((zt[2]-zt[1])/dense)*order+1
    local rix = {}
    local rix = {}
    local riy = {}
    local rox = {}
    local roy = {}

    -- Nodes in inner side
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        rix[i] = ri*cos(ang)
        riy[i] = ri*sin(ang)
        rox[i] = ro*cos(ang)
        roy[i] = ro*sin(ang)
    end

    -- Construct mesh
    for i = 1, nd do
        local xc = {rix[i],   riy[i],   zt[1],
                    rix[i],   riy[i],   zt[2],
                    rix[i+1], riy[i+1], zt[1],
                    rix[i+1], riy[i+1], zt[2],
                    rox[i],   roy[i],   zt[1],
                    rox[i],   roy[i],   zt[2],
                    rox[i+1], roy[i+1], zt[1],
                    rox[i+1], roy[i+1], zt[2]}
        mesh:add_block_shape(numnr, numns, numnz,etype, order, xc)
    end

end

function mesh_steady(ri,ro,zt,nd,angi)
   mesh_steady1(ri,ro,zt,nd,angi)
   mesh_steady1(ri,ro,zt,nd,angi+pi/4)
end

function mesh_transition1(ri,rt,ro,zt,nd,angi)

    local numnz  = ceil((zt[2]-zt[1])/dense)*order+1
    local angi= angi or 0
    local rix = {}
    local riy = {}
    local rtx = {}
    local rty = {}
    local rttx= {}
    local rtty= {}
    local rox = {}
    local roy = {}

    -- Nodes in inner side
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        rix[i] = ri*cos(ang)
        riy[i] = ri*sin(ang)
    end

    -- Nodes in outer side
    for i = 1, 3*nd+1 do
        local ang = pi/4/nd/3*(i-1) + angi
        rox[i] = ro*cos(ang)
        roy[i] = ro*sin(ang)
    end

    -- Nodes in transition region
    for i = 1, nd do
        local ang = pi/4/nd*(i-1) + angi
        local ang1= pi/4/nd/3
        local ang2= ang1*2
        rtx[2*i-1] = rt*cos(ang+ang1)
        rty[2*i-1] = rt*sin(ang+ang1)
        rtx[2*i  ] = rt*cos(ang+ang2)
        rty[2*i  ] = rt*sin(ang+ang2)
    end
    for i = 1, nd do
        local ang = pi/4/nd*(i-1) + angi
        local ang1= pi/4/nd/3
        local ang2= ang1*2
        rttx[3*i-2]= ri*cos(ang)
        rtty[3*i-2]= ri*sin(ang)
        rttx[3*i-1]= rt*cos(ang+ang1)
        rtty[3*i-1]= rt*sin(ang+ang1)
        rttx[3*i-0]= rt*cos(ang+ang2)
        rtty[3*i-0]= rt*sin(ang+ang2)
    end
    local ang = pi/4 + angi
    rttx[3*nd+1] = ri*cos(ang)
    rtty[3*nd+1] = ri*sin(ang)

    -- Construct innerside piece
    for i = 1, nd do
        local xc = {rix[i]    , riy[i],      zt[1],
                    rix[i]    , riy[i],      zt[2],
                    rix[i+1]  , riy[i+1],    zt[1],
                    rix[i+1]  , riy[i+1],    zt[2],
                    rtx[2*i-1], rty[2*i-1],  zt[1],
                    rtx[2*i-1], rty[2*i-1],  zt[2],
                    rtx[2*i]  , rty[2*i],    zt[1],
                    rtx[2*i]  , rty[2*i],    zt[2]}
        mesh:add_block_shape(order+1, order+1, numnz, etype, order, xc)
    end

    -- Construct outerside piece
    for i = 1, nd*3 do
        local xc = {rttx[i],   rtty[i],   zt[1],
                    rttx[i],   rtty[i],   zt[2],
                    rttx[i+1], rtty[i+1], zt[1],
                    rttx[i+1], rtty[i+1], zt[2],
                    rox[i],    roy[i],    zt[1],
                    rox[i],    roy[i],    zt[2],
                    rox[i+1],  roy[i+1],  zt[1],
                    rox[i+1],  roy[i+1],  zt[2]}
        mesh:add_block_shape(order+1, order+1, numnz, etype, order, xc)
    end

end

function mesh_transition(ri,rt,ro,zt,nd,angi)
   mesh_transition1(ri,rt,ro,zt,nd,angi)
   mesh_transition1(ri,rt,ro,zt,nd,angi+pi/4)
end

function mesh_post_interface1(rix,riy,ri,ro,zt,nd,angi)

    local angi   = angi or 0
    local numnr  = ceil((ro-ri)/dense)*order+1
    local numns  = order+1
    local numnz  = ceil((zt[2]-zt[1])/dense)*order+1
    local rox = {}
    local roy = {}

    -- Nodes in inner side
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        rox[i] = ro*cos(ang)
        roy[i] = ro*sin(ang)
    end

    -- Construct mesh
    for i = 1, nd do
        local xc = {rix[i],   riy[i],   zt[1],
                    rix[i],   riy[i],   zt[2],
                    rix[i+1], riy[i+1], zt[1],
                    rix[i+1], riy[i+1], zt[2],
                    rox[i],   roy[i],   zt[1],
                    rox[i],   roy[i],   zt[2],
                    rox[i+1], roy[i+1], zt[1],
                    rox[i+1], roy[i+1], zt[2]}
        mesh:add_block_shape(numnr, numns, numnz, etype, order, xc)
    end

end

function mesh_post_interface(xo,yo,ri,rt1,zt,nd,angi)

    local ax2 = {}
    local ay2 = {}

    -- Mesh portion1
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        ax2[i] = xo + ri*cos(ang)
        ay2[i] = yo + ri*sin(ang)
    end
    mesh_post_interface1(ax2,ay2,ri,rt1,zt,nd,angi)

    -- Mesh portion2
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + pi/4 + angi
        ax2[i] = xo + ri*cos(ang)
        ay2[i] = yo + ri*sin(ang)
    end
    mesh_post_interface1(ax2,ay2,ri,rt1,zt,nd,angi+pi/4)

end

function mesh_inner_disk(xo,yo,ri,zt,nd,angi)

    local numnz  = ceil((zt[2]-zt[1])/dense)*order+1
    local numnr = ceil(0.5*ri/dense)*order+1

    -- Obtain coordinates for inner portion
    local xc = { 0.0, 0.0, zt[1],
                 0.0, 0.0, zt[2],
        0.5*ri*cos(angi+pi/2), 0.5*ri*sin(angi+pi/2), zt[1],
        0.5*ri*cos(angi+pi/2), 0.5*ri*sin(angi+pi/2), zt[2],
        0.5*ri*cos(angi     ), 0.5*ri*sin(angi     ), zt[1],
        0.5*ri*cos(angi     ), 0.5*ri*sin(angi     ), zt[2],
        0.4*ri*cos(angi+pi/4)*sqrt(2), 0.4*ri*sin(angi+pi/4)*sqrt(2), zt[1],
        0.4*ri*cos(angi+pi/4)*sqrt(2), 0.4*ri*sin(angi+pi/4)*sqrt(2), zt[2]
    }

    for i = 1,8 do
        xc[3*i-2] = xc[3*i-2] + xo
        xc[3*i-1] = xc[3*i-1] + yo
    end 

    -- Mesh innerportion
    mesh:add_block_shape(order*nd+1, order*nd+1, numnz, etype, order, xc)

    -- Obtain coordinates for outerportion
    local ax1 = {}
    local ay1 = {}
    local ax2 = {}
    local ay2 = {}

    -- Mesh outer1 portion
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + angi
        ax1[i] = (nd+1 - i)/nd * xc[13] + (i-1)/nd * xc[19]
        ay1[i] = (nd+1 - i)/nd * xc[14] + (i-1)/nd * xc[20]
        ax2[i] = xo + ri*cos(ang)
        ay2[i] = yo + ri*sin(ang)
    end
    for i = 1, nd do
        local nodes = {ax1[i],   ay1[i],   zt[1],
                       ax1[i],   ay1[i],   zt[2],
                       ax1[i+1], ay1[i+1], zt[1],
                       ax1[i+1], ay1[i+1], zt[2],
                       ax2[i],   ay2[i],   zt[1],
                       ax2[i],   ay2[i],   zt[2],
                       ax2[i+1], ay2[i+1], zt[1],
                       ax2[i+1], ay2[i+1], zt[2]}
        mesh:add_block_shape(order*nd+1, numnr, numnz, etype, order, nodes)
    end

    -- Mesh outer2 portion
    for i = 1, nd+1 do
        local ang = pi/4/nd*(i-1) + pi/4 + angi
        ax1[i] = (nd+1 - i)/nd * xc[19] + (i-1)/nd * xc[7]
        ay1[i] = (nd+1 - i)/nd * xc[20] + (i-1)/nd * xc[8]
        ax2[i] = xo + ri*cos(ang)
        ay2[i] = yo + ri*sin(ang)
    end
    for i = 1, nd do
        local nodes = {ax1[i],   ay1[i],   zt[1],
                       ax1[i],   ay1[i],   zt[2],
                       ax1[i+1], ay1[i+1], zt[1],
                       ax1[i+1], ay1[i+1], zt[2],
                       ax2[i],   ay2[i],   zt[1],
                       ax2[i],   ay2[i],   zt[2],
                       ax2[i+1], ay2[i+1], zt[1],
                       ax2[i+1], ay2[i+1], zt[2]}
        mesh:add_block_shape(order*nd+1, numnr, numnz, etype, order, nodes)
    end

end

function mesh_pml(rbd, rpml)

    local ndivpml = ceil((rpml-rbd)/dense)*order+1
    local ndivb   = ceil( rbd      /dense)*order+1
    local nddiv   = nd*order+1
    mesh:blocks3dn( {   rbd, rpml     }, {ndivpml}, 
                    {     0,  rbd     }, {nddiv}, 
                    { -rpml, -rbd,  0 }, {ndivpml, ndivb}, etype, order)
    mesh:blocks3dn( {     0,  rbd     }, {nddiv}, 
                    {   rbd, rpml     }, {ndivpml}, 
                    { -rpml, -rbd,  0 }, {ndivpml, ndivb}, etype, order)
    mesh:blocks3dn( {   rbd, rpml     }, {ndivpml}, 
                    {   rbd, rpml     }, {ndivpml},
                    { -rpml, -rbd,  0 }, {ndivpml, ndivb}, etype, order)

end

function mesh_pmls(rbd, rpml)

    local ndivpml = ceil((rpml-rbd)/dense)*order+1
    local ndivb   = ceil( rbd      /dense)*order+1
    local nddiv   = nd*order+1

    mesh:blocks3dn( {   -rpml, -rbd, 0, rbd, rpml     }, {ndivpml, nddiv, nddiv, ndivpml}, 
                    {     rbd,  rpml                  }, {ndivpml}, 
                    {   -rpml, -rbd, 0                }, {ndivpml, ndivb}, etype, order)
    mesh:blocks3dn( {   -rpml, -rbd, 0, rbd, rpml     }, {ndivpml, nddiv, nddiv, ndivpml}, 
                    {   -rpml, -rbd                   }, {ndivpml}, 
                    {   -rpml, -rbd, 0                }, {ndivpml, ndivb}, etype, order)

    mesh:blocks3dn( { -rpml, -rbd      }, {ndivpml}, 
                    { -rbd ,    0,  rbd}, {nddiv, nddiv}, 
                    { -rpml, -rbd,  0  }, {ndivpml, ndivb}, etype, order)
    mesh:blocks3dn( {  rbd , rpml      }, {ndivpml}, 
                    { -rbd ,    0,  rbd}, {nddiv, nddiv}, 
                    { -rpml, -rbd,  0  }, {ndivpml, ndivb}, etype, order)

end

