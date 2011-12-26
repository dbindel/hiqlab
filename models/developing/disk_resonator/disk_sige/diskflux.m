clear;
qcloseall;

% -- Set parameters
params = [];
params.hdisk = 1.6e-6;
params.rdisk = 40e-6;
params.dense = 1e-6/2;
params.order = 3;

% -- Forcing frequency
wforce = 45.04 * 2e6*pi;

% -- Plotting parameters
plotopt = [];
plotopt.xscale = 1e-6;
plotopt.yscale = 1e-6;

% -- Load mesh from Lua input file
mesh  = Mesh_load('diskmesh.lua', params);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
drive_pat = Mesh_get_drive_f(mesh,'force_pattern');
harmonic_state(mesh,drive_pat,wforce);

% -- Plot flux
figure(1);
Mesh_make_harmonic(mesh, wforce);
p = Mesh_get_x(mesh);
E = Mesh_mean_power(mesh);
I = 1:9:size(p,2);
quiver(p(1,I), p(2,I), E(1,I), E(2,I)); 
