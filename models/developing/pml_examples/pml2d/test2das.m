% Driver for simple PML test: block sitting on a half space.
% Version 1: Use Lua for everything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test2das.m,v 1.1 2006/05/01 08:20:49 tkoyama Exp $

%try
 
  % -- Assemble and solve linear system

  w = 1;
  mesh  = Mesh_load('test2das.lua');
  Mesh_make_harmonic(mesh, w);
  [M,K] = Mesh_assemble_mk(mesh);
  F     = Mesh_assemble_R(mesh);
  u     = -(K-w^2*M)\F;
  Mesh_set_u(mesh, u);

  % -- Plot PML parameters

  figure(1);
  opt = [];
  opt.cfields = 'stretch_function';
  opt.ufields = 1;
  opt.ncfields = 2;
  plotfield2d(mesh, opt);

  % -- Plot forced response

  figure(2);
  opt = [];
  opt.ufields = [];
  opt.cfields = [1];
  opt.axequal = 1;
  plotcycle2d(mesh, 1, opt);
  Mesh_delete(mesh);

%catch
if 0

  % -- Clean up and pass on error
  if real(mesh) ~= 0
    Mesh_delete(mesh);
  end
  error(lasterr);

end
