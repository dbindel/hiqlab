-- Surrounding circuitry
xn1 = -160e-6
xn2 = -120e-6
xn3 =  -80e-6
yn1 =  -40e-6
yn2 =    0e-6
yn3 =   40e-6

xn  = { xn1, yn1,
        xn2, yn1,
        xn3, yn1,
        xn1, yn2,
        xn2, yn2,
        xn3, yn2,
        xn1, yn3,
        xn2, yn3,
        xn3, yn3}
numnp  = table.getn(xn)/2
idx    = {}
idx[1] = mesh:add_node(xn,numnp);
for i = 2,numnp do
  idx[i] = idx[i-1] + 1;
end

-- Add discrete circuit elements
resf   = mesh:Resistor(Rf)
ress   = mesh:Resistor(Rs)
wire   = mesh:VIsrc()
eno_vs =mesh:add_element({idx[ 7], idx[ 8]}, wire)
mesh:add_element({idx[ 8], idx[ 9]}, ress)
eno_vb = mesh:add_element({idx[ 5], idx[ 4]}, wire)
mesh:add_element({idx[ 5], idx[ 6]}, wire)
eno_vf = mesh:add_element({idx[ 2], idx[ 1]}, wire)
mesh:add_element({idx[ 2], idx[ 3]}, resf)

-- Add electrode element
eno_d,idg_d = add_electrode('electrode',idx[3],plateg_d,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_b,idg_b = add_electrode('electrode',idx[6],plateg_b,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_s,idg_s = add_electrode('electrode',idx[9],plateg_s,
                            get_dim_scale('Lt')/get_dim_scale('L'))

-- Nodal boundary conditions
function voltage_bc(x,y)
   if mesheq(x,xn1) and mesheq(y,yn1) then
                                  return '   u',   0; end
   if mesheq(x,xn1) and mesheq(y,yn2) then
                                  return '   u',   0; end
   if mesheq(x,xn1) and mesheq(y,yn3) then
                                  return '   u',   0; end
end

-- Element boundary conditions
function dc_voltage_source(elemnum)
  if elemnum == eno_vb then return 'f', Vb_dc; end
end
function ac_voltage_source(elemnum)
  if elemnum == eno_vf then return 'f', Vf_ac; end
end

-- Nodal sensing and forcing functions
function voltage_node2(x,y)
  if mesheq(x,xn2) and mesheq(y,yn1) then
     return '   u',  1; end;
end
function voltage_node9(x,y)
  if mesheq(x,xn3) and mesheq(y,yn3) then
     return '   u',  1; end;
end

-- Element sensing and forcing functions
function electrode_sense_V(elnum)
  if elnum == eno_s then return 'u',   1; end
end
function electrode_sense_Q(elnum)
  if elnum == eno_s then return ' u',   1; end
end
function electrode_drive_V(elnum)
  if elnum == eno_d then return 'u',   1; end
end
function electrode_drive_Q(elnum)
  if elnum == eno_d then return ' u',   1; end
end
function electrode_base_V(elnum)
  if elnum == eno_b then return 'u',   1; end
end
function electrode_drive_Q(elnum)
  if elnum == eno_b then return ' u',   1; end
end
