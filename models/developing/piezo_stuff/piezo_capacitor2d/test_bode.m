clear params;
clear H;
qcloseall;

% -- Meshfile
%meshfile = 'piezo_capacitor2d_f_mesh.lua';
meshfile = 'piezo_capacitor2d_wa_f_mesh.lua';

% -- Parameters
params.f0      = 0;
params.order   = 2;
params.dense   = 5;

opt.shift      = 180e6;
opt.use_matlab = 0;
opt.coupled    = 'electromech_elastic';


opt.plot_mesh  = 1;
opt.wr_min     = 0.85;
opt.wr_max     = 1.15;
opt.w_ndiv     = 200;
opt.lstyle     = 'r*-';
opt.visualQ    = 1;
[H,freq] = driver_piezo_capacitor2d_f_bode(meshfile,params,opt.shift,opt);
