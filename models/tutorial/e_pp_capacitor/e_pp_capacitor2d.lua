-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define mesh related global parameters
order   = order or 1       -- Order of element
dense   = dense or 10.0e-6  -- Approximate element size
meshtol = dense/100        -- Default mesh tolerance

-- Define element type
reps     = reps or 160
etype_es = mesh:Electrostatic2d(dim_scales.eps*reps)

-- Geometry parameters
Cw = 30e-6
Ch =  8e-6
Ct =  2e-6
dim_scales.Lt = Ct

-- Circuit parameters
Vf_dc= Vf_dc or 100  -- Forcing Plate DC
Vs_dc= Vs_dc or   0  -- Sensing Plate DC

-- Mesh center block
mesh:blocks2d({ -Cw/2, Cw/2},{0, Ch},etype_es,order,3*dense,dense)

-- Tie mesh together
mesh:tie()

-- Set boundary conditions
if use_case=='eppc2d_using_globals_sta.lua' then

  dofile 'eppc2d_using_globals_sta.lua'
  mesh:set_globals()
  mesh:set_globals_bc(global_voltage)

elseif use_case=='eppc2d_using_electrodes_sta.lua' then

  dofile 'eppc2d_using_electrodes_sta.lua'
  mesh:set_globals()
  mesh:set_bc{voltage_bc_t,voltage_bc_b}

else

  function vol_bc_t(x,y)
      if meshbetween(x, -Cw/2, Cw/2) and mesheq( y,Ch) then
                                     return 'u', Vf_dc;end
  end
  function vol_bc_b(x,y)
       if meshbetween(x, -Cw/2, Cw/2) and mesheq( y, 0) then
                                     return 'u',   0;end
  end
  mesh:set_bc{vol_bc_t,vol_bc_b}

end

-- Nodal sensing and forcing functions
function sense_top_Q(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y,Ch) then
                                  return 'f',  1;end
end
function sense_bot_Q(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y, 0) then
                                  return 'f',  1;end
end

sense_top_Q2 = make_nodal_func(sense_top_Q, 'f')
sense_bot_Q2 = make_nodal_func(sense_bot_Q, 'f')
