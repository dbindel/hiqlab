clear;
qcloseall;
% -- Parameters
param.reps     = 160;
param.dense    = 1e-6;
param.order    = 3;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('em_pp_capacitor2d.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt)

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
opt.mkc=1;
w  = 1.917995e+08*2*pi;
%drive_pat = Mesh_get_drive_elements_f(mesh,'ac_voltage_source');
drive_pat = Mesh_get_vector(mesh,'ac_voltage_source2');
forced_state(mesh,drive_pat,w,opt);

% -- Plot forced response
figure(3);
opt.deform  = 1e-4/Mesh_get_scale(mesh,'L');
opt.axequal = 1;
opt.cfields = [1,2,3];
plotcycle2d(mesh,opt.deform,opt);

% -- Clean up
Mesh_delete(mesh);
