clear;
qcloseall;
% -- Forcing frequency
w = 1;

% -- Load mesh from Lua input file
mesh  = Mesh_load('pmlaxis.lua');

% -- Plot the mesh
qfigure(1);
plotmesh(mesh);
%axis equal;

% -- Plot the stretch function
qfigure(2);
popt.cfields = 'stretch_function';
popt.ncfields= 2;
popt.axequal = 1;
popt.cbias   = 1;
plotfield2d(mesh,popt);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
F     =-Mesh_assemble_R(mesh);
forced_state(mesh,F,w);

% -- Plot forced response
qfigure(3);
opt.axequal = 1;
plotcycle2d(mesh,1,opt);

% -- Clean up
Mesh_delete(mesh);
