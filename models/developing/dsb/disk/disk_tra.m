clear;
qcloseall;
% -- Paramters
param.f0    =20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('diskmesh.lua',param);

% -- Form driving and sensing patterns
wc        = 2.889953e+08*2*pi;
drive_pat = Mesh_get_vector(mesh,'bode_drive_func2');
sense_pat = Mesh_get_vector(mesh,'bode_sense_func2');

% -- Compute transfer function
opt.wr_min    = 0.99;
opt.wr_max    = 1.01;
opt.w_ndiv    = 500;
opt.kmax      = 0;
opt.realbasis = 0;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Show bode plot
figure(2);
bopt.visualQ = 1;
plot_bode(freq,H,bopt);
