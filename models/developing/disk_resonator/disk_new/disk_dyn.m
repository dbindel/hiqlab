clear;
qcloseall;
% -- Parameters
param.f0    =20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('diskmesh.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute frequencies, Qs, and eigenvectors
w0      = 2.889953e+08*2*pi;
nev     = 1;
[V,w,Q] = mechmode(mesh, w0, nev);
[w_s,I] = sort((w-w0*ones(nev,1)));

% -- Show modes
figure(2);
opt.deform  = 1e-6/Mesh_get_scale(mesh,'L');
opt.axequal = 1;
opt.cfields = [1,2];
opt.cbias   = 1;
opt.animate = 1;
plot_mode(mesh,V(:,I),w(I),opt);
