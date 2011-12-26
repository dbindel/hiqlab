-- Add nodes
xn1= 0
yn1=-10e-6
yn2= 20e-6
xp = {xn1, yn1,
      xn1, yn2}
idx={}
idx[1] = mesh:add_node(xp,2);
idx[2] = idx[1] + 1

-- Define global shape function for plates
function plateg_d(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y,Ch) then
                                  return 1;end
end
function plateg_s(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y, 0) then
                                  return 1;end
end

-- Add Electrode element
eno_d,idg_d = add_electrode(mesh, 'electrode',idx[2],plateg_d,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_s,idg_s = add_electrode(mesh, 'electrode',idx[1],plateg_s,
                            get_dim_scale('Lt')/get_dim_scale('L'))

-- Nodal boundary conditions
function voltage_bc_t(x,y)
   if mesheq(x,xn1) and mesheq(y,yn2) then return ' u', Vf_dc; end
end
function voltage_bc_b(x,y)
   if mesheq(x,xn1) and mesheq(y,yn1) then return ' u',   0; end
end

-- Element sensing and forcing functions
electrode_top_V2 = global_indicator(idg_d)
electrode_bot_V2 = global_indicator(idg_s)
electrode_top_Q2 = branch_indicator(0,eno_d)
electrode_bot_Q2 = branch_indicator(0,eno_s)
