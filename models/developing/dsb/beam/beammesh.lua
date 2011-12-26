-- HiQLab
-- Copyright (c): Regents of the University of California

-- Include function definition file
require 'common.lua'
require 'beamshapes.lua'

-- Decide whether to reduce
use_full = true;

-- Define physical dimension of mesh
mesh  = Mesh:new(2)

-- Define approx size of elements and order
dense = 0.5e-6    -- Approximate element size (for block generator)
order = 2         -- Order of elements

-- Define size of beam
l = l or 10e-6    -- Beam length
w = 2e-6          -- Beam width

-- Define element type
mat   = mesh:PMLElastic2d_planestrain('silicon2')

-- Define mesh using block command
mesh:blocks2d( { 0, l }, { -w/2.0, w/2.0 }, mat )

if use_full then

  -- Define boundary conditions
  clamp_boundary(function(x,y) return x == 0 end, 'ux', 'uy')

else

  elt = impose_beam_x(0,l,-w/2,w/2)
  mesh:set_globals_bc(function(id)
    if id == elt[1] or id == elt[2] or id == elt[3] then
      return 'u', 0
    end
  end)

end

