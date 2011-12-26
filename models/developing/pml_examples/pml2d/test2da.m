% Driver for simple PML test: block sitting on a half space.
% Version 1: Use Lua for everything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test2da.m,v 1.1 2006/05/01 08:20:44 tkoyama Exp $

try
 
  % -- Assemble and solve linear system

  w = 1;
  mesh  = Mesh_load('test2da.lua');
  Mesh_make_harmonic(mesh, w);
  [M,K] = Mesh_assemble_mk(mesh);
  F     = Mesh_assemble_R(mesh);
  u     = -(K-w^2*M)\F;
  Mesh_set_u(mesh, u);

  % -- Plot forced response

  plotcycle2d(mesh);
  Mesh_delete(mesh);

catch

  % -- Clean up and pass on error
  if real(mesh) ~= 0
    Mesh_delete(mesh);
  end
  error(lasterr);

end
