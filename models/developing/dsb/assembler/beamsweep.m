% HiQLab
% Copyright (c): Regents of the University of California
% $Id: beamsweep.m,v 1.1 2006/05/01 03:47:43 tkoyama Exp $
clear;
qcloseall;

w = [];
l = linspace(10e-6, 25e-6, 16);
for k = 1:length(l)
  param.l = l(k);
  mesh  = Mesh_load('beammesh.lua', param);
  [M,K] = Mesh_assemble_mk(mesh);
  M = (M + M.')/2;
  K = (K + K.')/2;
  w(:,k)  = sqrt(eigs(K,M,4, 'sm'))/2/pi;
  Mesh_delete(mesh);
end
plot(l,real(w)); axis tight
