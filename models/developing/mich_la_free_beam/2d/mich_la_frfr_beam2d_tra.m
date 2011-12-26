clear;
qcloseall;
% -- Parameters
param.f0    = 20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('mich_la_frfr_beam2d.lua',param);

% -- Form driving and sensing patterns
ted_block_mesh(mesh);
wc        = 9.592613e+06*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'drive_pattern');
sense_pat = Mesh_get_sense_u(mesh,'sense_pattern');

% -- Compute transfer function
opt.wr_min = 0.95;
opt.wr_max = 1.05;
opt.w_ndiv = 300;
opt.mkc    = 1;
opt.kmax   = 0;
opt.viewmatrices = 0;
opt.structurep   = 0;
opt.realbasis    = 0;
[H,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Show bode plot
figure(2);
plot_bode(freq,H);
