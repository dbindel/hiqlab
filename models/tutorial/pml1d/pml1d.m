clear;
qcloseall;
% -- Pass parameters into Lua file
params.order = 3;
params.f0    =40;

% -- Forcing frequency
w     = 1;

% -- Load mesh from Lua input file
mesh  = Mesh_load('pml1d.lua', params);

% -- Display the mesh
figure(1)
plotmesh(mesh);

% -- Plot the stretch function
figure(2);
popt.cfields = 'stretch_function';
plotfield1d(mesh,popt);

% -- Construct mass, stiffness matrix, 
%    forcing vector and solve for 
%    harmonic displacements.
F     =-Mesh_assemble_R(mesh);
forced_state(mesh,F,w);

% -- Plot forced response
figure(3);
opt.deform = 1;
opt.axis   = [ 0, 12, -2, 2];
plotcycle1d(mesh, opt.deform,opt);

% -- Clean up
Mesh_delete(mesh);
