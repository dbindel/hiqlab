clear;
qcloseall;
% -- Parameters
param.use_case  = 'using_electrodes_dyn.lua';
param.order    = 1;
param.dense    = 1.0e-6;

param.Vf_dc    = 100;
param.fill_gap = 1;
param.reps     = 160;

param.f0       =   0;
f0             =  20;

w              = 140.00e6*2*pi;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('dielectric_drive.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt);
U = Mesh_get_u(mesh);

% -- Reload mesh with PML
Mesh_delete(mesh);
param.f0  = f0;
[mesh, L] = Mesh_load('dielectric_drive.lua',param);
Mesh_set_u(mesh,U);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
opt.mkc   = 1;
drive_pat = Mesh_get_drive_elements_f(mesh,'ac_voltage_source');
harmonic_state(mesh,drive_pat,w,opt);

% -- Plot forced response
figure(1);
opt.deform  = 1e-6;
opt.axequal = 1;
opt.cfields = [2];
plotcycle2d(mesh,opt.deform,opt);
