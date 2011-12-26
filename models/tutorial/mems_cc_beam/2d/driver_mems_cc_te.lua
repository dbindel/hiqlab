beamf = loadfile 'mems_cc_te.lua'

t0 = os.time()
order  = order  or 3
dense  = dense  or 0.5e-6
method = method or 'full'
beamf()
mesh:initialize()

w0 = 88*2e6*pi 
t1 = os.time()
if method == 'full' then
  status,dr,di = ted_compute_eigs(mesh, w0, 1)
else
  status,dr,di = ted_compute_eigsp(mesh, w0, 1)
end
t2 = os.time()

printf('Shift w0         : %g MHz', w0/2e6/pi)
printf('Damped w (real)  : %g MHz', dr[1]/2e6/pi)
printf('Damped w (imag)  : %g MHz', di[1]/2e6/pi)
printf('Q                : %g', dr[1]/(2*di[1]))
printf('Compute time     : %g', os.difftime(t2,t1))
printf('Init time        : %g', os.difftime(t1,t0))

mesh:delete()
