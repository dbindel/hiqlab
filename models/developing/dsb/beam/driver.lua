beamf = loadfile 'beammesh.lua'

w0  = 0;  -- Frequency estimate
nev = 5;  -- Number of eigs
ncv = 10; -- Size of space (~2 nev)
dr  = {}; -- Real parts of eigs
di  = {}; -- Imag parts of eigs

-- Load mesh and compute eigs
beamf()
mesh:initialize()
compute_eigs(mesh, w0, nev, ncv, dr, di);

-- Print eigs
for k = 1,5 do
  print(k, ':', dr[k]/2e6/pi, 'MHz')
end
