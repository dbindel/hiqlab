% -- Parameters
param.order = 3;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('ted_ring.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute frequencies, Qs, and eigenvectors
w0      = 705.50e6*2*pi;
nev     = 5;
[V,w,Q] = tedmode(mesh, w0, nev);
freq = w/2/pi;
fprintf('First fundamental[Hz]:%d\n',real(freq(1)));

% -- Put obtained displacements back into the mesh
show_mode = 1;
Mesh_set_u(mesh,V(:,show_mode));

% -- Plot mode shape
figure(2);
opt.deform   = 1e-6/Mesh_get_scale(mesh,'L');
popt.axequal = 1;
plotfield2d(mesh,opt);

% -- Plot animation of modeshape
figure(3);
opt.cfields = [3];
plotcycle2d(mesh,opt.deform,opt);
