% -- PML parameter
param.f0    = 0;
param.order = 3;
param.dense = 5e-6;

% -- Load mesh from Lua input file
mesh  = Mesh_load('ted_ring.lua',param);

% -- Plot the mesh
figure(1);
plotmesh(mesh);
axis equal;

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
opt.mkc=1;
w      = 618.50e6*2*pi;
%w      = 705.30e6*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'bode_force_function');
harmonic_state(mesh,drive_pat,w,opt);

% -- Plot forced response
figure(3);
opt.deform   = 1e-6;
opt.axequal = 1;
opt.cfields = [2];
plotcycle2d(mesh,opt.deform,opt);
