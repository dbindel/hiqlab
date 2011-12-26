beamf = loadfile 'mems_cc_te.lua'

-- Set parameters for analysis
order = 1
dense = 0.5e-6
etype = 'planestrain'
wafer = '111'
angle = math.pi/4
beamf()

-- Start analysis
mesh:initialize()
print('Numid:',mesh:get_numid())

-- Compute shift from preliminary analysis
status,dr,di = ted_compute_eigs_mech(mesh, 0.0, 1);
shift   = dr[1]
print('Shift:',shift/1e6/2/pi)

-- Computing actual frequencies
mesh:assign_ids()
status,dr,di = ted_compute_eigs(mesh, shift, 1);
freq = dr[1]/2/pi
Q    = sqrt(dr[1]^2 + di[1]^2)/2/di[1]
print('Freq :',freq/1e6)
print('Q    :',Q)
print('Re(w):',dr[1])
print('Im(w):',di[1])
