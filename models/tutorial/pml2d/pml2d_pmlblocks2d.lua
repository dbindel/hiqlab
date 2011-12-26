require 'common.lua'

order   = 3    -- Order of element
meshtol = 1e-6 -- Mesh tolerance
densex = 1     -- Approximate element size
densey = 1

xrad = 12  -- Radius of domain (x dimension)
yrad = 12  -- Radius of domain (y dimension)
dpml =  3  -- Thickness of PML (x dimension)

f0   = f0 or 40  -- Maximum stretch parameter

a = a or 1 -- Magnitude of vertical pull in forcing
b = b or 0 -- Magnitude of tilt in forcing

mesh = Mesh:new(2);
etype = mesh:PMLElastic2d_planestrain{ E = 10, nu = 0.3, rho = 2 }

-- Define pml block by a function
bc_pmlfunc, stretch_function = 
  pml_blocks2d({-1,1},{0,0},2,2*xrad,yrad,dpml,f0,etype,order,densex,densey)
mesh:tie()
etype:set_stretch(stretch_function)

-- Define boundary condition
function bc_function(x,y)
  if mesheq( y, 0) and meshleq( abs(x), 1) then return 'uu', 0, a+b*x; end
end
mesh:set_bc{bc_function,bc_pmlfunc}
