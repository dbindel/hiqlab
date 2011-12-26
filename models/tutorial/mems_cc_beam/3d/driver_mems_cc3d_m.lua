beamf = loadfile 'mems_cc3d_m.lua'

-- Set parameters for analysis
order = 3--3
dense = 3*1e-6
wafer = '100'
angle = 0 --math.pi/4
beamf()

-- Start analysis
mesh:initialize()
numid = mesh:get_numid()
print('Numid:',numid)

-- Compute shift from preliminary analysis
status,dr,di = compute_eigs(mesh, 0.0, 1);
shift   = dr[1]
--shift   = 76e6*2*pi
print('Shift:',shift/1e6/2/pi)
mesh:delete()
