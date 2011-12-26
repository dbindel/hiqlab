clear;
qcloseall;
% -- Parameters
param.order        = 1;
param.dense        = 1.0e-6;

param.fill_gap     =  1;
param.reps         =160;
param.Vb_dc        =100;

param.f0           =   0;
f0                 =  20;

wc         = 140.020e6*2*pi;
opt.wr_min = 0.9998;
opt.wr_max = 1.0002;
opt.w_ndiv = 100;
opt.kmax   = 0;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('de_drive_2port.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt);
U = Mesh_get_u(mesh);

% -- Reload mesh with PML
Mesh_delete(mesh);
param.f0  = f0;
[mesh, L] = Mesh_load('de_drive_2port.lua',param);
Mesh_set_u(mesh,U);

% -- Compute frequencies, Qs, and eigenvectors
drive_pat = Mesh_get_drive_elements_f(mesh,'ac_voltage_source');
sense_pat1= Mesh_get_sense_u(mesh,'voltage_node2');
sense_pat2= Mesh_get_sense_u(mesh,'voltage_node9');
sense_pat = [sense_pat1,sense_pat2];

opt.mkc  =  1;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Extract variables
Rf = Lua_get_double(L,'Rf');
Rs = Lua_get_double(L,'Rs');
H = H(2,:)./H(1,:)*(Rs+Rf)/Rs;

% -- Show bode plot
figure(1);
opt.lstyle    = 'r';
opt.magnitude = 1;
plot_bode(freq,H,opt);
