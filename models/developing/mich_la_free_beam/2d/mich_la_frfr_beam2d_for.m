clear;
qcloseall;
% -- PML parameter
param.f0    = 20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
mesh  = Mesh_load('mich_la_frfr_beam2d.lua',param);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
opt = [];
opt.mkc   = 1;
w         = 9.588685e+06*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'drive_pattern');
forced_state(mesh,drive_pat,w,opt);

% -- Plot forced response
figure(3);
opt.deform  = 1e-7;
opt.axequal = 1;
opt.cfields = [1,2,3];
opt.cbias   =  1;
plotcycle2d(mesh,opt.deform,opt);
