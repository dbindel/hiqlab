beamf = loadfile 'mems_cc3d_te.lua'

-- Set parameters for analysis
order = 3
dense = 3*1e-6
wafer = '100'
angle = 0 --math.pi/4
beamf()

-- Start analysis
mesh:initialize()
numid = mesh:get_numid()
print('Numid:',numid)

-- Compute shift from preliminary analysis
status,dr,di = ted_compute_eigs_mech(mesh, 0.0, 1)
shift   = dr[1]
print('Shift:',shift/1e6/2/pi)

-- Computing actual frequencies
mesh:assign_ids()
status,dr,di = ted_compute_eigs(mesh, shift, 1)
freq = dr[1]/2/pi
Q    = sqrt(dr[1]^2 + di[1]^2)/2/di[1]
print('Shift:',shift/1e6/2/pi)
print('Freq :',freq/1e6)
print('Q    :',Q)
print('Re(w):',dr[1])
print('Im(w):',di[1])
print('Numid:',numid)
mesh:delete()
