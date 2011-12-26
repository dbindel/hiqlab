-- Define mesh function
function connect_beam(quad, r_mid, ang_l, ang_r, bl_l, bl_r, bw_l, bw_r, etype, order, dense)

    local xl1,xl2,xr1,xr2
    xr1= r_mid*cos(ang_r)
    xr2= xr1 + bl_r
    xl1=-r_mid*cos(ang_l)
    xl2= xl1 - bl_l

    if quad==1 then
      mesh:blocks2d({xr1,xr2},{ 0,  bw_r/2},etype,order,dense)
    elseif quad==4 then
      mesh:blocks2d({xr1,xr2},{ -bw_r/2, 0},etype,order,dense)
    elseif quad==2 then
      mesh:blocks2d({xl2,xl1},{  0, bw_l/2},etype,order,dense)
    elseif quad==3 then
      mesh:blocks2d({xl2,xl1},{-bw_l/2,  0},etype,order,dense)
    end 
end

function connect_piece(quad, r_in, r_mid, ang_l, ang_r, ndiv_t, bw_l, bw_r, etype, order, dense)

    local rat_l = (r_mid*sin(ang_l))/(bw_l/2)
    local rat_r = (r_mid*sin(ang_r))/(bw_r/2)

    local ndiv_bl1= order*ceil(          bw_l/2/dense)
    local ndiv_bl2= order*ceil((rat_l-1)*bw_l/2/dense)
    local ndiv_br1= order*ceil(          bw_r/2/dense)
    local ndiv_br2= order*ceil((rat_r-1)*bw_r/2/dense)
    local div1,div2

    local xl1,xl2,xr1,xr2
    local yl1,yl2,yr1,yr2

    xl1=-r_in *cos(ang_l)
    yl1= r_in *sin(ang_l)
    xl2=-r_mid*cos(ang_l)
    yl2= r_mid*sin(ang_l)
    xr1= r_in *cos(ang_r)
    yr1= r_in *sin(ang_r)
    xr2= r_mid*cos(ang_r)
    yr2= r_mid*sin(ang_r)

    local nodex1, nodex2
    if quad==1 or quad==4 then
      nodex1 = { xr1,         0, xr1,  yr1/rat_r,   xr2,         0, xr2, yr2/rat_r}
      nodex2 = { xr1, yr1/rat_r, xr1,  yr1      ,   xr2, yr2/rat_r, xr2, yr2      }
      if quad==4 then
          nodex1 = reflect_quad(nodex1,{0,1})
          nodex2 = reflect_quad(nodex2,{0,1})
      end
      div1 = ndiv_br1
      div2 = ndiv_br2
    elseif quad==2 or quad==3 then
      nodex1 = { xl2,         0, xl2, yl2/rat_l,    xl1,         0, xl1, yl1/rat_r}
      nodex2 = { xl2, yl2/rat_l, xl2, yl2      ,    xl1, yl1/rat_r, xl1, yl1      }
      if quad==3 then
          nodex1 = reflect_quad(nodex1,{0,1})
          nodex2 = reflect_quad(nodex2,{0,1})
      end
      div1 = ndiv_bl1
      div2 = ndiv_bl2
    end

    mesh:add_block_shape(ndiv_t+1, div1+1, etype, order, nodex1)
    mesh:add_block_shape(ndiv_t+1, div2+1, etype, order, nodex2)

end


function layer_ring2d(quad, r_in, r_out, ang_l, ang_r, ndiv_r, ndiv_t, etype, order)

    local ang1, ang2

    if quad==1 then 
        ang1 = pi/2
        ang2 = ang_r
    elseif quad==2 then
        ang1 = pi-ang_l
        ang2 = pi/2
    elseif quad==3 then
        ang1 = pi*3/2
        ang2 = pi+ang_l
    elseif quad==4 then
        ang1 =2*pi-ang_r
        ang2 =  pi*3/2
    end
    function ring(x,y)
        local anew = ang1 + (x+1)/2 * (ang2 -ang1);
        local rnew = r_in + (y+1)/2 * (r_out-r_in);
        return rnew*cos(anew),rnew*sin(anew);
    end
    mesh:add_block_transform(ndiv_r+1, ndiv_t+1, etype, order, ring)
end

function quad_ring2d(quad, r_in, r_t, slit_t, slit_l, slit_r, bl_l, bl_r, bw_l, bw_r, etype, order, dense)

    local r_out = r_in  + r_t
    local r_mid = r_out - slit_t
    local ang_l = asin(slit_l/2/r_mid)
    local ang_r = asin(slit_r/2/r_mid)

    local ndiv_ti= order*ceil((r_mid-r_in)/dense)  -- Along straight inner edge
    local ndiv_to= order*ceil( slit_t     /dense)  -- Along straight outer edge
    local ndiv_r = 10*ndiv_ti                      -- Along quater of ring perimeter
       
    layer_ring2d(quad,r_in, r_mid, ang_l, ang_r, 
                 ndiv_r, ndiv_ti, etype, order)
    if slit_t~=0 then
    layer_ring2d(quad,r_mid, r_out, ang_l, ang_r, 
                 ndiv_r, ndiv_to, etype, order)
    end
    
    connect_piece(quad, r_in, r_mid, ang_l, ang_r, ndiv_ti, bw_l, bw_r, etype, order, dense)
--    connect_beam(quad, r_mid, ang_l, ang_r, bl_l, bl_r, bw_l, bw_r, etype, order, dense)

end

    

