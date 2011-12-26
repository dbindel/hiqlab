clear;
qcloseall;

% -- Meshfile
meshfile = 'piezobeam2d_mesh.lua';
%meshfile = 'piezobeam2d_wa_mesh.lua';

% -- Parameters
params.f0      = 20;
params.order   = 3;
params.dense   = 1;

opt.numeigs    = 1;
opt.eigs_n     = 7;
opt.shift      = 184e6;
opt.use_matlab = 0;
opt.coupled    = 'piezoelectric_elastic';

opt.kmax       = 10;

opt.plot_defo  = 1;
opt.plot_mesh  = 1;
opt.plot_mag   = 1e-1;
opt.cfields    = [1];

% Test
%driver_piezoelectric_elastic_mode(meshfile,params,opt);
%driver_piezoelectric_elastic_force(meshfile,params,opt.shift,opt);
[H,freq] = driver_piezoelectric_elastic_bode(meshfile,params,opt.shift,opt);
