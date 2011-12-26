-- Include function definition file
require 'common.lua'
require 'quad_ring2d.lua'

-- Define physical dimension of mesh
mesh = Mesh:new(2)

-- Define nondimensionalization parameters
mtype = 'poly_silicon'

-- Define element type
etype = mesh:PMLElastic2d_planestress(mtype)

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

-- Functions
rot_ang = rot_ang or 0
function transform(x,y)
    local x1,y1
    x1 = cos(rot_ang)*x - sin(rot_ang)*y
    y1 = sin(rot_ang)*x + cos(rot_ang)*y
    return x1,y1
end
function inv_transform(x,y)
    local x1, y1
    x1 =  cos(-rot_ang) * (x) - sin(-rot_ang) * (y)
    y1 =  sin(-rot_ang) * (x) + cos(-rot_ang) * (y)
    return x1,y1
end
function reflect_y(x,y)
    return x,-y
end
function reflect_x(x,y)
    return -x,y
end
function reflect_conn(ix_old)
    local ix_new = {}
    for i = 0,order do
        for j = 0,order do
             ix_new[i*(order+1)+1+j] = ix_old[(order-i)*(order+1)+1+j]
--             ix_new[i*(order+1)+1+j] = ix_old[i*(order+1)+1+order-j]
        end
    end
    return ix_new
end
-- Set boundary conditions
function bc_func(x1,y1)
--    x,y = x1,y1
    x,y = inv_transform(x1,y1)
--    if mesheq(x,-Bl-Rm*cos(Al)) or mesheq(x, Bl+Rm*cos(Al)) then
--                                return 'uu', 0, 0; end
    if (mesheq(x,-Rm*cos(Al)) or mesheq(x, Rm*cos(Al))) 
                             and meshbetween(y,-Bw/2,Bw/2) then
                                return 'uu', 0, 0; end
    local rad = sqrt(x^2+y^2)
    if mesheq(rad,Ri+Rt) and meshgeq(abs(y),(Ri+Rt)*sin(Al)) then
       local xf,yf = transform(x/rad,y/rad)
       return 'ff', xf, yf;
--       return 'ff', x/rad, y/rad;
    end
end

-- Construct ring
nstart = mesh:numnp()
quad_ring2d(3, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(4, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(1, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
quad_ring2d(2, Ri, Rt, Nt, Nw, Nw, Bl, Bl, Bw, Bw, etype, order, dense)
nend = mesh:numnp()
mesh:transform_nodes(nstart,nend-1,transform)

-- Set boundary conditions
mesh:set_bc{bc_func}

-- Tie mesh together
mesh:tie()

-- Set forcing and sensing function
function bode_force_function(x1,y1)
    x,y = inv_transform(x1,y1)
    local rad = sqrt(x^2+y^2)
    if mesheq(rad,Ri+Rt) and meshgeq(abs(y),(Ri+Rt)*sin(Al)) then
       local xf,yf = transform(x/rad,y/rad)
       return 'ff', xf, yf;
--       return 'ff', x/rad, y/rad;
    end
end
function bode_sense_function(x,y)
   local rad = sqrt(x^2+y^2)
   if mesheq(rad,Ri+Rt) and meshgeq(abs(y),(Ri+Rt)*sin(Al)) then
      return 'uu', x/rad, y/rad;
   end
end

