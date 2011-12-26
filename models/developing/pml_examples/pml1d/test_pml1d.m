% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_pml1d.m,v 1.1 2006/05/01 08:20:31 tkoyama Exp $

mesh  = Mesh_load('pml1d.lua');
k = 2*pi/10;
Mesh_make_harmonic(mesh, k);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
u     = -(K-k^2*M)\F;
Mesh_set_u(mesh, u);

x = Mesh_get_x(mesh);
u = Mesh_get_disp(mesh);
plot(x, real(u), x, imag(u));

Mesh_delete(mesh);
