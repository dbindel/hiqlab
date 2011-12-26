% HiQLab
% Copyright (c): Regents of the University of California
% $Id: beamsweep.m,v 1.1 2006/05/01 03:47:43 tkoyama Exp $
clear;
close all;

nu = linspace(0, 0.48, 50);
for k = 1:length(nu)
  param.nu = nu(k);
  mesh  = Mesh_load('beammesh2.lua', param);
  [M,K] = Mesh_assemble_mk(mesh);
  M = (M + M.')/2;
  K = (K + K.')/2;
  w(:,k)  = sqrt(eigs(K,M,14, -0.1))/2/pi;
  Mesh_delete(mesh);
end
plot(nu,real(w(1:end-3,:)), '.'); axis tight
