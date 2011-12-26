% Scan deformed shape in an interesting frequency range

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_bode2.m,v 1.2 2006/05/04 05:37:49 tkoyama Exp $

smin =  90e6 * 2*pi;
smax = 100e6 * 2*pi;
s0   = 100e6 * 2*pi;
k    = 100;
n    = 100;

mesh  = Mesh_load('checkermesh.lua');
[M,K] = Mesh_assemble_mk(mesh);
L     = Mesh_get_sense_u(mesh, 'sense_checker');
B     = Mesh_get_drive_f(mesh, 'drive_checker');

freq = linspace(smin, smax, n);
[Mk,Kk,Lk,Bk,Vk] = rom_arnoldi(M,K,L,B, k, s0);
H = compute_bode_mech(freq, Mk,Kk,Lk,Bk);

for i = 1:n
 
  uk = (Kk - freq(i)^2 * Mk)\Bk;
  u  = Vk * uk;
  u  = u * 4e-6 / max(abs(u));

  Mesh_set_u(mesh, u);
  opt.deform = 1;
  plotmesh(mesh,opt);

  refresh;
  pause(0.05);

end

Mesh_delete(mesh);
