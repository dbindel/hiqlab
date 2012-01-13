% Driver for simple PML test: block sitting on a half space.
% Version 1: Use Lua for everything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_marray.m,v 1.4 2006/05/02 04:15:21 tkoyama Exp $
clear;
  
% --- Assemble and solve linear system

w = 1;
[mesh,L] = Mesh_load('pml2d.lua');
Mesh_make_harmonic(mesh, w);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
A     = K-w^2*M;
u     = -(A\F);

% --- Compare to the same computation done in Lua
Lua_set_array(L, 'A', A);
Lua_set_arrayz(L, 'F', F);
Lua_dofile(L, 'test_marray.lua');

AA = Lua_get_array(L, 'A');
FF = Lua_get_array(L, 'F');
x  = Lua_get_array(L, 'x');

norm(x + u)

Mesh_delete(mesh);

