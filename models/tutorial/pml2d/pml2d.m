clear;
qcloseall;
% -- Forcing frequency
w = 1;

% -- Forcing pattern
params.a = 1;
params.b = 0;

% -- PML parameter
params.f0 =40;

% -- Load mesh from Lua input file
mesh  = Mesh_load('pml2d.lua',params);
%mesh  = Mesh_load('pml2d_pmlblocks2d.lua',params);

% -- Plot the mesh
qfigure(1);
plotmesh(mesh);

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
opt.subplot = [1,2];
plotcycle2d(mesh,1,opt);

% -- Clean up
Mesh_delete(mesh);
