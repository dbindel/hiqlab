-- Define global shape function for plates
function top_g(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y,Ch) then
                                  return 1;end
end
function bot_g(x,y)
   if meshbetween(x, -Cw/2, Cw/2) and mesheq( y, 0) then
                                  return 1;end
end

-- Add global variables
idg_t = mesh:add_global(top_g,'V','Q')
idg_b = mesh:add_global(bot_g,'V','Q')

-- Define global variable boundary condition function
function global_voltage(idg)
  if idg == idg_t then return 'u', Vf_dc; end
  if idg == idg_b then return 'u',   0; end
end
