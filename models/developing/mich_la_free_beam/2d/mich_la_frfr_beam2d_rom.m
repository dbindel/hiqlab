clear;
qcloseall;
% -- Parameters
param.f0    = 20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('mich_la_frfr_beam2d.lua',param);

% -- Compute frequencies, Qs, and eigenvectors
ted_block_mesh(mesh);
wc        = 9.592613e+06*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'force_pattern');
sense_pat = Mesh_get_sense_u(mesh,'sense_pattern');

opt.wr_min = 0.95;
opt.wr_max = 1.05;
opt.w_ndiv = 100;
opt.mkc    = 1;
opt.viewmatrices = 0;

% -- Full case
opt.kmax   = 0;
opt.structurep = 0;
opt.realbasis  = 0;
[H,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Reduced model
opt.kmax   = 8;
opt.structurep = 0;
opt.realbasis  = 0;
[H1,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Real basis
opt.kmax   = 4;
opt.structurep = 0;
opt.realbasis  = 1;
[H2,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Structure basis
opt.kmax   = 4;
opt.structurep = 1;
opt.realbasis  = 0;
[H3,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Real and structure basis
opt.kmax   = 2;
opt.structurep = 1;
opt.realbasis  = 1;
[H4,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

clf;
figure(1);
plot(freq/1e6/2/pi,log10(abs(H1-H)./abs(H)),'r*-');
hold on;
plot(freq/1e6/2/pi,log10(abs(H2-H)./abs(H)),'b*-');
plot(freq/1e6/2/pi,log10(abs(H3-H)./abs(H)),'g*-');
plot(freq/1e6/2/pi,log10(abs(H4-H)./abs(H)),'k*-');
xlabel('Frequency [MHz]');
ylabel('Relative error');
