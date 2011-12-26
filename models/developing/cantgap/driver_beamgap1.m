% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_beamgap1.m,v 1.1 2005/06/28 19:05:44 dbindel Exp $

function [Vall, fall] = driver_beamgap1;

print_newton = 1;
popt.deform = 1;
popt.cfields = 3;

[mesh,Ls] = Mesh_load('beamgap.lua'); 
numid = Mesh_get_numid(mesh);

u  = zeros(numid,1);

Vall = 10:60;
lall = zeros(1,length(Vall));

for iV=1:length(Vall)

  % Set new voltage
  VV = Vall(iV);
  fprintf('Voltage %e\n', VV);
  Lua_set_double(Ls, 'V', VV);
  Mesh_apply_bc(mesh);

  % Modified Newton corrector
  K = Mesh_assemble_dR(mesh, 1,0,0);
  for k = 1:8 
    R  = Mesh_assemble_R(mesh);
    du = K\R;
    u  = u-du;
    Mesh_set_u(mesh, u);
    if print_newton
      fprintf('  %d: %e %e\n', k, norm(du), norm(R));
    end
  end

  % Show the mesh
  plotfield2d(mesh, popt);
  shg;

  % Compute eigenvalues on reduced system
  [M,K] = Mesh_assemble_mk(mesh);
  [L,U,P,Q] = qlu(K);
  opts.disp = 0;
  d = eigs(@Kfun,numid, 3, 'LM', opts, L,U,P,Q, M);
  f = sqrt(1/d(1)) / 2 / pi;
  fall(iV) = f;
  fprintf('  Freq: %e kHz\n', f);

end

Mesh_delete(mesh);

% Plot voltage vs frequency
figure(2);
plot(Vall, fall);
xlabel('Potential (V)');
ylabel('Frequency (kHz)');


% -- Apply K\M

function [y] = Kfun(x, L, U, P, Q, M)

y = Q*(U\(L\(P*(M*x))));
