-- Surrounding circuitry
xn1 =  -60e-6
xn2 =   60e-6
xn3 =  100e-6
xn4 =  140e-6
xn5 =  260e-6

yn1 =  -80e-6
yn2 =  -40e-6
yn3 =    0e-6

xn  = { xn1, yn1,
        xn1, yn2,
        xn1, yn3,
        xn2, yn3,
        xn3, yn1,
        xn3, yn2,
        xn3, yn3,
        xn4, yn3,
        xn5, yn1,
        xn5, yn3}
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

eno_vf = mesh:add_element({idx[ 2], idx[ 1]}, wire)
mesh:add_element({idx[ 2], idx[ 3]}, resf)

eno_vb = mesh:add_element({idx[ 6], idx[ 5]}, wire)
mesh:add_element({idx[ 6], idx[ 7]}, ress)

eno_vs =mesh:add_element({idx[10], idx[ 9]}, wire)

mesh:add_element({idx[ 4], idx[ 7]}, wire)
mesh:add_element({idx[ 7], idx[ 8]}, wire)

-- Add electrode element
eno_d1,idg_d1 = add_electrode('electrode',idx[3],plateg_d1,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_s1,idg_s1 = add_electrode('electrode',idx[4],plateg_s1,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_d2,idg_d2 = add_electrode('electrode',idx[8],plateg_d2,
                            get_dim_scale('Lt')/get_dim_scale('L'))
eno_s2,idg_s2 = add_electrode('electrode',idx[10],plateg_s2,
                            get_dim_scale('Lt')/get_dim_scale('L'))

-- Nodal boundary conditions
function voltage_bc(x,y)
   if mesheq(x,xn1) and mesheq(y,yn1) then
                                  return '   u',   0; end
   if mesheq(x,xn3) and mesheq(y,yn1) then
                                  return '   u',   0; end
   if mesheq(x,xn5) and mesheq(y,yn1) then
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
  if mesheq(x,xn1) and mesheq(y,yn2) then
     return '   u',  1; end;
end
function voltage_node7(x,y)
  if mesheq(x,xn3) and mesheq(y,yn3) then
     return '   u',  1; end;
end
