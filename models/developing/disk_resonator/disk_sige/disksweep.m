% Driver to compare simulated frequencies and Q values with
% measured values reported in the Michigan Transducers 03 paper
% on diamond disk.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: disksweep.m,v 1.1 2006/05/01 04:59:41 tkoyama Exp $

param = []; 
param.order = 3;
param.dense = 1e-6/3;
w0 = 2*pi*47.5e6;
t1 = linspace(1.505e-6,  1.53e-6, 10);
t2 = linspace(1.53e-6,   1.55e-6, 20);
t3 = linspace(1.55e-6,   1.6e-6,  10);
t = [t1(1:end-1), t2(1:end-1), t3];
wall = [];
Qall = [];

for k = 1:length(t)
  param.hdisk = t(k);
  mesh = Mesh_load('diskmesh.lua', param);
  [M,K] = Mesh_assemble_mk(mesh);
  Mesh_delete(mesh);

  fprintf('hdisk = %d: %d/%d\n', t(k), k, length(t));
  [VV,w,Q] = pml_mode(M,K,w0,2);
  [Q2k,I] = sort(Q);
  Qall(:,k) = Q(I);
  wall(:,k) = w(I) / 2e6/pi;
end

semilogy(t/1e-6, Qall); axis tight
xlabel('Film thickness (um)');
ylabel('Q');
