-- Surrounding circuitry
xn1 =  -75e-6
xn2 =   75e-6
yn1 = -100e-6
yn2 =  -25e-6
yn3 =   0.0

xn  = { xn1, yn1,
        xn1, yn2,
        xn1, yn3,
        xn2, yn1,
        xn2, yn2,
        xn2, yn3 }
numnp  = table.getn(xn)/2
idx    = {}
idx[1] = mesh:add_node(xn,numnp);
for i = 2,numnp do
  idx[i] = idx[i-1] + 1;
end

-- Add electrode element
eno_d,idg_d = add_electrode('electrode',idx[3],plateg_d,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_s,idg_s = add_electrode('electrode',idx[6],plateg_s,
                            get_dim_scale('Lt')/get_dim_scale('L'))

-- Nodal boundary conditions
function voltage_bc(x,y)
   if mesheq(x,xn1) and mesheq(y,yn3) then
                                  return '   u', Vf_dc; end
   if mesheq(x,xn2) and mesheq(y,yn3) then
                                  return '   u',   0; end
end

-- Element sensing and forcing functions
function electrode_top_V(elnum)
  if elnum == eno_d then return 'u',   1; end
end
function electrode_top_Q(elnum)
  if elnum == eno_d then return ' u',   1; end
end
function electrode_bot_V(elnum)
  if elnum == eno_s then return 'u',   1; end
end
function electrode_bot_Q(elnum)
  if elnum == eno_s then return ' u',   1; end
end
