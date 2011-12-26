-- Michigan Free-Free Beam Example
-- Article:
-- Wan-Thai Hsu, John R. Clark, and Clark T.-C. Nguyen
-- Q-Optimized Lateral Free-Free Beam Micromechanical Resonators
-- Digest of Technical Papers, the 11th Int. Conf. on Solid State
--  Sensors and Actuators(Transducers'01), Munich Germany
-- pp.1110--1113

-- Include function definition file
require 'common.lua'
require 'generate_frfr_beam2d.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
ted_nondim('silicon2',2e-6)

-- Define mesh related global parameters
order   = order or 1       -- Order of element
dense   = dense or 0.5e-6  -- Approximate element size
meshtol = dense/100        -- Default mesh tolerance

-- Define element type
mtype = 'poly_silicon'
etype = mesh:PMLElastic2d_te_planestress(mtype)

-- Define geometry of domain
Bl = Bl or 40.4e-6 -- Beam length
Bw = Bw or  2e-6   -- Beam width
Bt = Bt or  2e-6   -- Beam thickness
Sl = Sl or 25.6e-6 -- Support beam length
Sw = Sw or 1.2e-6  -- Support beam width
Sx = Sx or 10e-6   -- Support beam location(from center)
Dw = Dw or  4e-6   -- Drive width

Pw = Pw or 12e-6   -- Anchor width
Ph = Ph or  6e-6   -- Anchor height
Pd = Pd or  3e-6   -- PML depth    
f0 = f0 or  20

-- Generate Free-Free beam
bc_func,st_func = generate_frfr_beam2d( 0e-6, 0e-6,   0)

-- Set stretch function
function st_funcs(x,y)
    local sx = 0
	local sy = 0
	for i=1, 4 do
	    local tx, ty = st_func[i](x,y)
		sx = sx + tx
		sy = sy + ty
	end
	return sx, sy
end

etype:set_stretch(st_funcs)

-- Tie mesh together
mesh:tie()

-- Set boundary conditions
mesh:set_bc(bc_func)

--Define driving and sensing functions
function drive_pattern(x,y)
   -- Forcing pattern
   if (mesheq(y,-Bw/2) or mesheq(y, Bw/2)) and meshbetween(x,-Dw/2,Dw/2) then
                                           return ' f ', 1/Bt;
   end
end
function sense_pattern(x,y)
   -- Sensing pattern
   if (mesheq(y,-Bw/2) or mesheq(y, Bw/2)) and meshbetween(x,-Dw/2,Dw/2) then
                                           return ' u ', 1;
   end
end
