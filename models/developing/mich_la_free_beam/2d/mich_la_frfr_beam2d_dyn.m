clear;
qcloseall;
% -- Parameters
param.f0    =  0;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('mich_la_frfr_beam2d.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute frequencies, Qs, and eigenvectors
w0      = 9.6055e6*2*pi;
nev     = 1;
[V,w,Q] = tedmode(mesh, w0, nev);

% -- Show modes
figure(2);
opt.deform  = 1e-4/Mesh_get_scale(mesh,'L');
opt.axequal = 1;
opt.cfields = [3];
opt.cbias   = 1;
opt.animate = 1;
plot_mode(mesh,V,w,opt);
