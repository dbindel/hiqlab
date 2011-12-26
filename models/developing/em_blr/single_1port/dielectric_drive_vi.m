clear;
qcloseall;
% -- Parameters
param.use_case  = 'using_electrodes_vi.lua';

param.order        = 1;
param.dense        = 1.0e-6;

param.fill_gap     = 1;
param.reps         = 160;
param.Vf_dc        = 100;

param.f0           =   0;
f0                 =  20;

wc             = 140.020e6*2*pi;
opt.wr_min     = 0.9995;
opt.wr_max     = 1.0005;
opt.w_ndiv     = 100;
opt.kmax       =  0;

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

% -- Compute transfer function
drive_pat = Mesh_get_drive_elements_f(mesh,'ac_voltage_source');
sense_pat = Mesh_get_sense_elements_u(mesh,'current_voltage_source');
opt.mkc   =  1;
[H,freq]  = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
H         =-H/Lua_get_double(L,'Vf_ac');

% -- Show bode plot
figure(2);
opt.lstyle    = 'r';
opt.magnitude = 1;
plot_bode(freq,H,opt);
