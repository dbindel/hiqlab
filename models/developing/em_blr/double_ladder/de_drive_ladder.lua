-- Include function definition file
require 'common.lua'
require '../single_1port/blr.lua'

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
cap_mat  = get_material(cap_mat)
gap_mat  = 'hfo2'
fill_gap = fill_gap or 1
reps     = reps     or 160
cap_mat.rho = cap_mat.rho * 1.0
etype_c1 = mesh:PMLELastic2d_planestress(cap_mat)
cap_mat.rho = cap_mat.rho * 1.00002
etype_c2 = mesh:PMLElastic2d_planestress(cap_mat)
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
Vf_dc= Vf_dc or   0  -- Forcing Plate DC
Vb_dc= Vb_dc or 100  -- Forcing Plate DC
Vs_dc= Vs_dc or   0  -- Sensing Plate DC
Vf_ac= Vf_ac or 0.1  -- Forcing Plate AC
Vb_ac= Vb_ac or   0  -- Forcing Plate AC
Vs_ac= Vs_ac or   0  -- Sensing Plate AC
Rf   = Rf    or 50  -- Resistor(Ohms)force
Rs   = Rs    or 50  -- Resistor(Ohms)sense
Rf   = 1000
Rs   = 10000

-- Construct BLR1
xc    = 0
yc    = 0
alpha = 0
port  = 1
etype_c = etype_c1
bc_func1,st_func1,gs_func1 = construct_blr(port,xc,yc,alpha)
plateg_d1 = gs_func1[1]
plateg_s1 = gs_func1[2]

-- Construct BLR2
xc    =200e-6
yc    = 0
alpha = 0
port  = 1
etype_c = etype_c2
bc_func2,st_func2,gs_func2 = construct_blr(port,xc,yc,alpha)
plateg_d2 = gs_func2[1]
plateg_s2 = gs_func2[2]

-- Use electrode elements
require 'using_electrodes_ladder.lua'

-- Construct tables for boundary conditions and 
--                      stretch functions
bc_func = {}
st_func = {}
for i = 1,table.getn(bc_func1) do
    bc_func[i] = bc_func1[i]
    st_func[i] = st_func1[i]
end
num_bc = table.getn(bc_func)
for i = 1,table.getn(bc_func2) do
    bc_func[num_bc+i] = bc_func2[i]
    st_func[num_bc+i] = st_func2[i]
end
bc_func[table.getn(bc_func)+1] = voltage_bc


-- Set stretch function
etype_c1:set_stretch(st_func)
etype_c2:set_stretch(st_func)

-- Tie mesh together
mesh:tie()

-- Set global variables and related boundary conditions
mesh:set_globals()
mesh:set_elements_bc(dc_voltage_source)

-- Set boundary conditions
mesh:set_bc(bc_func)
