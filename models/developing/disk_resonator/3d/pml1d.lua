-- Include function definition file
require 'common.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(1)

-- Define element type
etype = mesh:PMLScalar1d(1, 1);

-- Define order of element
order = order or 1

-- Define mesh tolerance
meshtol = 1e-6

-- Define geometry of domain 
xrad     =  4              -- Radius of domain (x dimension)
xpml     =  3              -- Start of PML (x dimension)
num_elem = 12              -- Elements in x
--xrad     = 12              -- Radius of domain (x dimension)
--xpml     =  9              -- Start of PML (x dimension)
--num_elem = 36              -- Elements in x

-- PML parameter
f0   = f0 or  40  -- Maximum stretch parameter

-- Define mesh using block command
mesh:add_block( 0, xrad, num_elem+1, etype, order)

-- Define stretching function
function stretch_function(x)

  local sx = 0
  sx = max( (abs(x)-xpml)/(xrad-xpml) * f0, 0)
  
  return sx
end

-- Set stretch function to element
etype:set_stretch(stretch_function)

-- Define boundary condition
function bc_function(x)
  if mesheq( x,    0) then return 'u', 1 end
  if mesheq( x, xrad) then return 'u', 0 end
end

mesh:set_bc(bc_function)
