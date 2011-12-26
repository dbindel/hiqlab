-- HiQLab
-- Copyright (c): Regents of the University of California

 
beamf = loadfile 'beamgap.lua'
beamf();
mesh:initialize();

numid = mesh:get_numid();
u  = QArray:new(numid, 1);
R  = QArray:new(numid, 1);
du = QArray:new(numid, 1);

-- FIXME: I need a way to move the no-mass directions up front.

for VV=10,60 do

  -- Set new voltage
  printf('Voltage %e', VV);
  V = VV;
  mesh:apply_bc()

  -- Modified Newton corrector
  K = mesh:assemble_dR(1,0,0);
  for k = 1,8 do
    mesh:assemble_R();
    mesh:get_f(R);
    K:solve(du,R);
    u:sub(du);
    mesh:set_u(u)
    if print_newton then
      printf('  %d: %e %e', k, du:normf(), R:normf());
    end
  end
  K:delete()

  -- Compute eigenvalues on reduced system
  K = mesh:assemble_dR(1,0,0);
  M = mesh:assemble_dR(0,0,1);

  d,v = leigs{
    op = function(x,y)
      scratch = M:apply(x,scratch)
      K:solve(y,scratch)
    end,
    n = numid,
    real = true,
    mode = 'standard',
    which = 'LM',
    nev = 1
  }
  printf('  Freq: %e kHz', sqrt(1/d(1)) / 2 / pi);

  K:delete()
  M:delete()
end
