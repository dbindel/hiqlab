clear;
qcloseall;

% -- Meshfile
meshfile = 'piezoring2hd_mesh.lua';

% -- Parameters
params.order = 3;
params.dense = 10;  % Does not make difference in this mesh
                    % Only for meshtol

opt.numeigs  = 10;
opt.whichmode=  1;
opt.coupled  = 'elastic';
opt.kmax     = 20;
%opt.shift    = 232.0e6;
opt.shift    = 720.5e6;
wforce       = opt.shift;

opt.plot_defo= 10;
opt.plot_mesh= 1;
opt.plot_mag = 1e-3;
opt.cfields  = [1];
opt.nframes  = 16;
opt.wr_min   = 0.95;
opt.wr_max   = 1.05;
opt.w_ndiv   = 200;

% Test 
driver_piezoelectric_elastic_mode(meshfile,params,opt);
%driver_piezoelectric_elastic_force(meshfile,params,wforce,opt);
%driver_piezoelectric_elastic_bode(meshfile,params,wforce,opt);
