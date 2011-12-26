-- Include function definition file
require 'common.lua'
require 'quad_ring2d.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
mtype = 'poly_silicon'
--mtype = 'sc_silicon'
ted_nondim(mtype,1e-6)

-- Define element type
caxis    = { 1, 0, 0, 0, 1, 0}           -- Crystal orientation
--etype = make_material_te(mtype, 'planestress',caxis)
etype = mesh:PMLELastic2d_te_planestress(mtype)

-- Define mesh related global parameters
order   = order or 1         -- Order of element
dense   = dense or 10.0e-6   -- Approximate element size
meshtol = dense/100          -- Default mesh tolerance

-- Define geometry of domain
Nw = Nw or 10e-6   -- Notch width
Nt = Nt or 10e-6   -- Notch thickness
Ri = Ri or 90e-6   -- Ring inner radius
Rt = Rt or 20e-6   -- Ring thickness
Rm = Ri+Rt-Nt      -- Ring middle radius
Bw = Bw or 6e-6    -- Beam width
Bl = Bl or 15e-6   -- Beam length
Pd = Pd or 10e-6   -- PML depth
Pw = Pw or 50e-6   -- Anchor width
Ph = Ph or 25e-6   -- Anchor height
f0 = f0 or 20
Al = asin(Nw/2/Rm) -- Angle of the notch opening

-- Construct ring
quad_ring2d(1, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(2, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(3, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(4, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)

-- Construct pmlblocks
bc_pml1,s_pml1 =  pml_blocks2d('E',{-Bl-Rm*cos(Al)},
                                   {-Bw/2,0,Bw/2},3,Pw,Ph,Pd,f0,etype,order,dense)
bc_pml2,s_pml2 =  pml_blocks2d('W',{ Rm*cos(Al)+Bl },
                                   {-Bw/2,0,Bw/2},3,Pw,Ph,Pd,f0,etype,order,dense)

-- Set stretch function
etype:set_stretch{s_pml1,s_pml2}

-- Set boundary conditions
mesh:set_bc{bc_pml1,bc_pml2}

-- Tie mesh together
mesh:tie()

-- Set forcing and sensing function
function bode_force_function(x,y)
   local rad = sqrt(x^2+y^2)
   if mesheq(rad,Ri+Rt) and meshgeq(abs(y),(Ri+Rt)*sin(Al)) then
      return 'ff ', x/rad, y/rad;
   end
end
function bode_sense_function(x,y)
   local rad = sqrt(x^2+y^2)
   if mesheq(rad,Ri+Rt) and meshgeq(abs(y),(Ri+Rt)*sin(Al)) then
      return 'uu', x/rad, y/rad;
   end
end

