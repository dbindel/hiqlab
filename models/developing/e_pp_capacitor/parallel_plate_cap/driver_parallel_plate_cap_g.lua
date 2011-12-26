-- Parallel plate capacitance test calculation
-- Uses Global shape functions
 
-- HiQLab
-- Copyright (c): Regents of the University of California

capf = loadfile 'parallel_plate_cap_g.lua'
capf()
mesh:initialize()

-- Set the total charge (dual to the tied voltage) to 1
numid = mesh:get_numid()
xglb  = mesh:globalid(idg1)+1
qq = QArray:new(numid,1)
qq:set(xglb,1, 1)

-- Find the corresponding voltage distribution
K  = mesh:assemble_dR(1,0,0)
vv = K:solve(qq)

-- Compute the capacitance as Q/V and by the hand formula
C1 = eps0*w/g
C2 = qq(xglb)/vv(xglb)

print('Capacitance by hand:   ', C1)
print('Capacitance from mesh: ', C2)

mesh:delete()
