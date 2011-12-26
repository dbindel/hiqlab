-- Add global variables
idg1 = mesh:add_global(plateg_d,'V','Q')
idg2 = mesh:add_global(plateg_s,'V','Q')

-- Set global boundary conditions
function global_voltage(idg)
  if idg == idg1 then return 'u', Vf_dc; end
  if idg == idg2 then return 'u',     0; end
end
