-- Parallel plate capacitance test calculation
-- Uses TieField element
 
-- HiQLab
-- Copyright (c): Regents of the University of California

capf = loadfile 'parallel_plate_cap_f.lua'
capf()
mesh:initialize()

-- Set the total charge (dual to the tied voltage) to 1
numid = mesh:get_numid()
xtie  = mesh:branchid(0,ntie)+1
qq = QArray:new(numid,1)
qq:set(xtie,1, 1)

-- Find the corresponding voltage distribution
K  = mesh:assemble_dR(1,0,0)
vv = K:solve(qq)

-- Compute the capacitance as Q/V and by the hand formula
C1 = eps0*w/g
C2 = qq(xtie)/vv(xtie)

print('Capacitance by hand:   ', C1)
print('Capacitance from mesh: ', C2)

mesh:delete()
