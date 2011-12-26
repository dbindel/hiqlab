-- HiQLab
-- Copyright (c): Regents of the University of California

-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh  = Mesh:new(2)

-- Define approx size of elements and order
dense = 0.5e-6    -- Approximate element size (for block generator)
order = 2         -- Order of elements

-- Define size of beam
l = 10e-6         -- Beam length
w = 2e-6          -- Beam width

-- Define element type
mat   = mesh:PMLElastic2d_planestrain('silicon2')

-- Define mesh using block command
mesh:blocks2d( { 0, l }, { -w/2.0, w/2.0 }, mat )

-- Define boundary conditions
clamp_boundary(function(x,y) return x == 0 end, 'ux', 'uy')
