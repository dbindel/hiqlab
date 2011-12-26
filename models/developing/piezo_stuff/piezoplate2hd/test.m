clear;
qcloseall;

% Test 1
params.order   = 3;
params.dense   = 0.05;
params.numeigs = 1;
params.plotdefo= 1;
params.plotmesh= 1;
%params.shift   = 23.3e6;
params.shift   = 82.0e6;
%params.shift   = 97.5e6;
wforce   = params.shift;
%driver_piezoelectric_elastic('piezoplate2hd_blr_mesh.lua',params);
%pause;

% Test 2
qcloseall;
opt.cfields = [1];
opt.plot_mag=1e-2;
opt.nframes = 16;
%wforce   = 82.0e6;
%wforce   = 97.5e6;
driver_piezoelectric_elastic_force('piezoplate2hd_blr_mesh.lua',params,wforce,opt);
pause;

% Test 3
qcloseall;
wforce   = 23.3e6;
%wforce   = 82.0e6;
%wforce   = 97.5e6;
driver_piezo_bode('piezoplate2hd_blr_mesh.lua',params,wforce);
