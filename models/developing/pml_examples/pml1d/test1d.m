% Driver for simple PML test: block sitting on a half space.
% Version 1: Use Lua for everything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test1d.m,v 1.1 2006/05/01 08:20:34 tkoyama Exp $

% -- Assemble and solve linear system

param.k     = 5;
param.order = 3;
param.w     = 1;

mesh  = Mesh_load('test1d.lua', param);
Mesh_make_harmonic(mesh, param.w);
[M,K] = Mesh_assemble_mk(mesh);
F     =-Mesh_assemble_R(mesh);
u     = (K-param.w^2*M)\F;
Mesh_set_u(mesh, u);

% -- Plot forced response

opt.axis = [0, 12, -1, 1];
plotcycle1d(mesh, 1, opt);
Mesh_delete(mesh);
