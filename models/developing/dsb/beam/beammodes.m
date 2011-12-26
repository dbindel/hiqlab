% HiQLab
% Copyright (c): Regents of the University of California
% $Id: beammodes.m,v 1.1 2006/05/01 03:47:43 tkoyama Exp $

clear;
qcloseall;
% -- Load mesh from Lua input file
mesh = Mesh_load('beammesh.lua');

% -- Display the mesh
figure(1);
opt.axequal = 1; 
plotmesh(mesh,opt); 

% -- Assemble matrices and solve eigenvalue problem
[M,K] = Mesh_assemble_mk(mesh);
M = (M + M.')/2;
K = (K + K.')/2;

if size(M) > 10
  [V,D] = eigs(K,M, 4, 'sm');
  w = sqrt(diag(D))/2/pi;
  [w,I] = sort(w);
  V     = V(:,I);
  show_mode = 1;
else
  [V,D] = eig(full(K),full(M));
  w = sqrt(diag(D))/2/pi;
  [w,I] = sort(w);
  V     = V(:,I);
  show_mode = 1;
end

% -- Put obtained displacements back into the mesh
Mesh_set_u(mesh,V(:,show_mode));

% -- Display displacements
figure(2);
opt.axequal = 1;
opt.deform = 1e-10;
plotfield2d(mesh, opt);

% -- Delete mesh
Mesh_delete(mesh);
