% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_taxis.m,v 1.1 2005/03/11 20:01:16 dbindel Exp $

param.ltype = 0;
nmodes = 10;

mesh = Mesh_load('test_taxis.lua', param);

if param.ltype == -1
  maxdof = 2;
else
  maxdof = 3;
end

figure(1);
plotmesh(mesh); axis equal

[M,K] = Mesh_assemble_mk(mesh);
M = (M + M.')/2;
K = (K + K.')/2;
[V,D] = eigs(K,M, nmodes, -0.1);
w = sqrt(diag(D))/2/pi;
[w,I] = sort(real(w));
V = V(:,I);

for k = 1:nmodes
  figure(2);
  opt.axequal = 1;
  opt.deform = 1;
  opt.ufields = [1 maxdof]; 
  opt.cfields = 2;
  opt.cbar = 1;
  opt.titles = sprintf('omega = %g', w(k));
  Mesh_scale_u(mesh, V(:,k), 1:maxdof, 1);
  plotfield2d(mesh, opt);
  pause;
end

Mesh_delete(mesh);
