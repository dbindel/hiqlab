beamf = loadfile 'mich_la_frfr_beam2d.lua'

-- Set parameters for analysis
f0   =  0.0
order= 3
dense= 1.0e-6

-- Frequency extraction parameters
w0  = 9.6055e6*2*pi;
nev = 1

-- Load mesh and compute eigs
beamf()
mesh:initialize()
print('Numid:',mesh:get_numid())
status,dr,di,vr,vi = ted_compute_eigs(mesh, w0, nev);

---[[
mesh:set_u(vr)
mesh:set_ui(vi)
dxf = DXFile:new('mich2d')
dxf:writemesh(mesh)
dxf:delete()
--]]

-- Print eigs
for k = 1,nev do
  Q    = sqrt(dr[k]^2 + di[k]^2)/2/di[k]
  print(k, ':', dr[k]/2e6/pi, '  ', di[k]/2e6/pi, 'MHz', '  Q:', Q)
end
  
mesh:delete()
