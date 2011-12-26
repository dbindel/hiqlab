% Scan deformed shape in an interesting frequency range

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_gui.m,v 1.2 2006/05/04 05:37:52 tkoyama Exp $

smin =  90e6 * 2*pi;
smax =  98e6 * 2*pi;
s0   = 100e6 * 2*pi;
k    = 100;
n    = 200;

mesh  = Mesh_load('checkermesh.lua');
[M,K] = Mesh_assemble_mk(mesh);
L     = Mesh_get_sense_u(mesh, 'sense_checker');
B     = Mesh_get_drive_f(mesh, 'drive_checker');

freq = linspace(smin, smax, n);
[Mk,Kk,Lk,Bk,Vk] = rom_arnoldi(M,K,L,B, k, s0);
H = compute_bode_mech(freq, Mk,Kk,Lk,Bk);

driver_gui(0, mesh, Mk, Kk, Bk, Vk, freq, H);

