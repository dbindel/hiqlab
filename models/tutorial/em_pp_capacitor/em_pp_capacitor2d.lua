-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
em_nondim('silicon2',2e-6)

-- Define mesh related global parameters
order   = order or 1       -- Order of element
dense   = dense or 10.0e-6  -- Approximate element size
meshtol = dense/100        -- Default mesh tolerance

-- Define element type
cap_mat  = 'silicon2'
reps     = reps or 160
etype_c  = mesh:PMLElastic2d_planestress(cap_mat)
etype_em = mesh:CoupleEM2d(dim_scales.eps*reps)

-- Geometry parameters
Cw = 30e-6
Ch =  8e-6
Ct =  2e-6
dim_scales.Lt = Ct

-- Circuit parameters
Vf_dc= Vf_dc or 100  -- Forcing Plate DC
Vs_dc= Vs_dc or   0  -- Sensing Plate DC
Vf_ac= Vf_ac or 0.1  -- Forcing Plate AC
Vs_ac= Vs_ac or   0  -- Sensing Plate AC

-- Mesh center block
mesh:blocks2d({ -Cw/2, Cw/2},{-Ch/2, Ch/2},etype_c ,order,3*dense,dense)
mesh:blocks2d({ -Cw/2, Cw/2},{-Ch/2, Ch/2},etype_em,order,3*dense,dense)

-- Define mechanical boundary condition
clamp_boundary( 
  function(x,y) return meshbetween(x, -Cw/2,Cw/2) and mesheq(y,-Ch/2) end,
  'ux', 'uy')

-- Tie mesh together
mesh:tie()

-- Set boundary conditions
dofile 'emppc2d_using_electrodes.lua'
mesh:set_globals()
