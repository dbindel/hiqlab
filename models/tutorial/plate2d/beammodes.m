% HiQLab
% Copyright (c): Regents of the University of California
% $Id: beammodes.m,v 1.1 2006/05/01 03:47:43 tkoyama Exp $

%clear;
%close all;
% -- Load mesh from Lua input file
param.nu = 0.2;
mesh = Mesh_load('beammesh2.lua', param);

% -- Display the mesh
figure(1);
opt.axequal = 1; 
%plotmesh(mesh,opt); 

% -- Assemble matrices and solve eigenvalue problem
[M,K] = Mesh_assemble_mk(mesh);
M = (M + M.')/2;
K = (K + K.')/2;
[V,D] = eigs(K,M, 10, 'sm');
w = sqrt(diag(D))/2/pi

[w,I] = sort(w);
V     = V(:,I);

for show_mode = 4:11

% -- Put obtained displacements back into the mesh
%show_mode = 6;
Mesh_set_u(mesh,V(:,show_mode));

% -- Display displacements
figure(show_mode);
opt.axequal = 1;
opt.deform = 1;
plotfield2d(mesh, opt);

end

% -- Delete mesh
Mesh_delete(mesh);
