-- Include function definition file
require 'common.lua'
require 'blr.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
em_nondim('silicon2',2e-6)

-- Define mesh related global parameters
order   = order or 1       -- Order of element
dense   = dense or 2.0e-6  -- Approximate element size
meshtol = dense/100        -- Default mesh tolerance

-- Define element type
cap_mat  = 'silicon2'
gap_mat  = 'hfo2'
fill_gap = fill_gap or 1
reps     = reps     or 160

etype_c  = mesh:PMLELastic2d_planestress(cap_mat)
etype_g  = mesh:PMLELastic2d_planestress(gap_mat)
etype_em = mesh:CoupleEM2d(dim_scales.eps*reps)

-- Geometry parameters
-- Beam
Bh = 90.4e-6
Bw = 12e-6
Lt =  2e-6
dim_scales.Lt = Lt
-- Connection beam
Cw =  0.5e-6
Cl =  8e-6
--Csh= 'straight'
Csh= 'serpentine'
-- Dielectric material
Dt = 0.2e-6
Dx =  30e-6
-- PML
Pw =  30e-6
Ph =  15e-6
Pd =   5e-6
f0 =  f0 or 20

-- Circuit parameters
Vf_dc= Vf_dc or 100  -- Forcing Plate DC
Vs_dc= Vs_dc or   0  -- Sensing Plate DC
Vf_ac= Vf_ac or 0.1  -- Forcing Plate AC
Vs_ac= Vs_ac or   0  -- Sensing Plate AC
Rf   = Rf    or 50   -- Resistor(Ohms)force side
Rs   = Rs    or 50   -- Resistor(Ohms)sense side

-- Construct BLR
xc    = 0
yc    = 0
alpha = 0
port  = 1
bc_func,st_func,gs_func = construct_blr(port,xc,yc,alpha)
plateg_d = gs_func[1]
plateg_s = gs_func[2]

if use_case == 'using_globals_sta.lua' then

    -- Use global variables
    require 'using_globals_sta.lua'
    mesh:set_globals()
    mesh:set_globals_bc(global_voltage)

elseif use_case == 'using_electrodes_sta.lua' then

    -- Use electrode elements
    require 'using_electrodes_sta.lua'
    bc_func[table.getn(bc_func)+1] = voltage_bc
    mesh:set_globals()

elseif use_case == 'using_electrodes_dyn.lua' then

    -- Use electrode elements
    require 'using_electrodes_dyn.lua'
    bc_func[table.getn(bc_func)+1] = voltage_bc
    mesh:set_globals()
    mesh:set_elements_bc(dc_voltage_source)

elseif use_case == 'using_electrodes_vi.lua' then

    -- Use electrode elements
    require 'using_electrodes_vi.lua'
    bc_func[table.getn(bc_func)+1] = voltage_bc
    mesh:set_globals()
    mesh:set_elements_bc(dc_voltage_source)

else 
    error('No support')
end

-- Set stretch function
etype_c:set_stretch(st_func)

-- Tie mesh together
mesh:tie()

-- Set boundary conditions
mesh:set_bc(bc_func)

-- Set nodal sensing and forcing functions
function sense_top_Q(x,y)
   if meshbetween(x, -Bw/2, Bw/2) and (mesheq( y,Dx+Dt) or mesheq(y,-Dx-Dt)) then
                                  return '  f',  1;end
end
function sense_bot_Q(x,y)
   if meshbetween(x, -Bw/2, Bw/2) and (mesheq( y,Dx) or mesheq(y,-Dx)) then
                                  return '  f',  1;end
end
function sense_gap_dist(x,y)
   if mesheq(x, -Bw/2) and mesheq( y,Dx+Dt)  then
                                  return ' u ', -1;end
   if mesheq(x, -Bw/2) and mesheq( y,Dx)  then
                                  return ' u ',  1;end
end
