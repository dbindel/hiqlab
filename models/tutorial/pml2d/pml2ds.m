clear;
qcloseall;
% -- Forcing frequency
w = 3;

% -- Load mesh from Lua input file
mesh  = Mesh_load('pml2ds.lua');

% -- Plot the mesh
qfigure(1);
plotmesh(mesh);
% axis equal;

% -- Plot the stretch function
qfigure(2);
popt.cfields = 'stretch_function';
popt.ufields = 2;
popt.ncfields= 2;
popt.axequal = 1;
popt.cbias   = 1;
plotfield2d(mesh,popt);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
%    U = K\F;
F     =-Mesh_assemble_R(mesh);
forced_state(mesh,F,w);

% -- Plot forced response
qfigure(3);
opt.axequal = 1;
opt.ufields = 1;
opt.cfields = 1;
plotcycle2d(mesh,0,opt);

% -- Clean up
Mesh_delete(mesh);
