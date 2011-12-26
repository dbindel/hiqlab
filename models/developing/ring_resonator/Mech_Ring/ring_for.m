clear
qcloseall;
% -- PML parameter
param.f0    = 0;
param.order = 3;
param.dense = 5e-6;
param.rot_ang = 0;

% -- Load mesh from Lua input file
%mesh  = Mesh_load('ring_only.lua',param);

%mesh  = Mesh_load('ring_only.lua',param);

%param.which_half = 'W';
%mesh  = Mesh_load('ring_half_only.lua',param);

%param.which_quad = 1;
%mesh  = Mesh_load('ring_quad_only.lua',param);

mesh  = Mesh_load('bar.lua',param);

% -- Plot the mesh
figure(1);
plotmesh(mesh);
axis equal;

% -- Quiver plot of forces
figure(2);
node_force = Mesh_get_sense_force(mesh,'bode_force_function');
x          = Mesh_get_x(mesh);
quiver(x(1,:),x(2,:),node_force(1,:),node_force(2,:));
axis equal;

pause;

% -- Static displacement
figure(3);
static_state(mesh);
popt.axequal = 1;
plotfield2d(mesh,popt);
pause;

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
opt.mkc=1;
w      = 618.50e6*2*pi;
%w      = 705.30e6*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'bode_force_function');
harmonic_state(mesh,drive_pat,w,opt);

% -- Plot forced response
figure(4);
opt.deform   = 1e-6;
opt.axequal = 1;
opt.cfields = [1,2,3];
plotcycle2d(mesh,opt.deform,opt);
