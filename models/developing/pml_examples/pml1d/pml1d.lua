require 'common.lua'

epw   = epw   or 10  -- Elements per wave
dpw   = dpw   or 1   -- Damping per wave
Ne1   = Ne1   or 10  -- Elements in first region
Ne2   = Ne2   or 10  -- Elements in second region
order = order or 3   -- Polynomial order

-- Assemble mesh quantities
mesh = Mesh:new(1);
D    = mesh:PMLScalar1d(1,1);
mesh:add_block(-Ne1, Ne2, order*(Ne1+Ne2)+1, D, order);

-- Define stretch and BC
D:set_stretch(function(x)
   return max(x*dpw/epw,0)
end)

mesh:set_bc(function(x)
   if x == -Ne1 then return 'u', 1 end  
end)
