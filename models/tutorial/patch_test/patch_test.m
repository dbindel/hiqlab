clear;
qcloseall;
% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('patch_test.lua');

% -- Display the mesh
qfigure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Construct stiffness matrix and forcing vector
%    and solve for displacements
%    U = K\F;
static_state(mesh);

% -- Display displacements
qfigure(2);
opt.deform = 10;
plotfield2d(mesh,opt);

% -- Clean up
Mesh_delete(mesh);
