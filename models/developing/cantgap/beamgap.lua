-- HiQLab
-- Copyright (c): Regents of the University of California

require 'common.lua'

l = 100           -- Beam length
g = 2             -- Gap width
w = 2             -- Beam width
dense = 0.25       -- Approximate element size (for block generator)
order = 1         -- Order of elements
eps0  = 8.854e-6  -- Permittivity of free space
V     = V or 20   -- Voltage applied

E   = 165e3;
nu  = 0.3;
rho = 2300e-12;

mesh   = Mesh:new(2)
mat    = mesh:PMLElastic2d(E, nu, rho, 0);
mat_em = mesh:CoupleEM2d(eps0);
mesh:blocks2d( { 0, l }, { g, g+w }, mat )
mesh:blocks2d( { 0, l }, { 0, g }, mat_em, order, dense , dense/2 )
mesh:tie(1e-3)

function anchor_bc(x,y)
--  if (x == 0) and (y >= g) then return 'uu', 0, 0; end
  if (y == g+w)            then return 'uu', 0, 0; end
  if (y == 0)              then return 'uu', 0, 0; end
end

function gap_x_bc(x,y) if (y <  g) then return  'uu',0,  0;     end end
--function gap_v_bc(x,y) if (y <= g) then return '  u', V*y/g; end end
function gap_x1_bc(x,y) if (y == 0) then return '  u',      0;  end end
function gap_x2_bc(x,y) if (y == g) then return '  u',      V;  end end

mesh:set_bc{anchor_bc, gap_x_bc, gap_x1_bc, gap_x2_bc}
--mesh:set_bc{anchor_bc, gap_x_bc, gap_v_bc}
