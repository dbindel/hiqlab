beamf = loadfile 'ffbeam3d_mesh.lua'

-- Set parameters for analysis
order = 2
wafer = '100'
angle = 0
rdense= 10--0.2
rf0   = 20.0

-- Start preliminary analysis(No PML)
--f0      = 0.0
--dense   = 10
--beamf()
--mesh:initialize()

--status,dr,di  = ted_compute_eigs_mech(mesh, 0.0, 1);
shift         = 9.6e6*2*pi
print('Shift:',shift/1e6/2/pi)
--mesh:delete()


-- Computing actual frequencies using shift
f0    = rf0
dense = rdense
beamf()
mesh:initialize()
print('Numid:',mesh:get_numid())

status,dr,di = ted_compute_eigs(mesh, shift, 1);
freq = dr[1]/2/pi
Q    = sqrt(dr[1]^2 + di[1]^2)/2/di[1]
print('Freq :',freq/1e6)
print('Q    :',Q)
print('Re(w):',dr[1])
print('Im(w):',di[1])
mesh:delete()
