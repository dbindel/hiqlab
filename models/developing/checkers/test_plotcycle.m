% Do a movie of the motion excited at 95 MHz

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_plotcycle.m,v 1.3 2006/05/04 05:37:51 tkoyama Exp $

s0   = 95e6 * 2*pi;

mesh  = Mesh_load('checkermesh.lua');
[M,K] = Mesh_assemble_mk(mesh);
L     = Mesh_get_sense_u(mesh, 'sense_checker');
B     = Mesh_get_drive_f(mesh, 'drive_checker');

u = (K-s0^2*M)\B;
u = 3e-6 * u / max(abs(u));
Mesh_set_u(mesh, u);

opt = [];
opt.cfields = [];
plotcycle2d(mesh,1,opt);
Mesh_delete(mesh);
