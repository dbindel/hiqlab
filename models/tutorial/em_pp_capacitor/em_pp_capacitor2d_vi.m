clear;
qcloseall;
% -- Parameters
param.dense = 1e-6;
param.reps  = 160;
param.order = 3;

wc  = 1.917998e+08*2*pi;
opt.wr_min = 0.99999;
opt.wr_max = 1.00001;
opt.w_ndiv = 200;
opt.kmax   = 0;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('em_pp_capacitor2d.lua',param);

% -- Compute for static state
sopt.nonlinear = 'NR';
static_state(mesh,sopt)

% -- Form driving and sensing pattern
drive_pat = Mesh_get_vector(mesh,'ac_voltage_source2');
sense_pat = Mesh_get_vector(mesh,'sense_motional_current2');
Vac       = Lua_get_double(L,'Vf_ac');

% -- Compute transfer function
opt.mkc   = 1;
[H,freq]  = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
H = H(1,:)/Vac;

% -- Show bode plot
figure(2);
opt.lstyle = 'r*-';
plot_bode(freq,H,opt);

% -- Clean up
Mesh_delete(mesh);
