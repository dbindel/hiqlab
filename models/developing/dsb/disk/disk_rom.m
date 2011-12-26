%clear;
%qcloseall;
% -- Paramters
param = [];
param.f0    =20;
param.order = 3;
param.dense = 1e-6;
%param.beam_assumption = 1;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('diskmesh.lua',param);

% -- Compute frequencies, Qs, and eigenvectors
wc        = 2.889953e+08*2*pi;
drive_pat = Mesh_get_vector(mesh,'bode_drive_func2');
sense_pat = Mesh_get_vector(mesh,'bode_sense_func2');

opt.wr_min    = 0.99;
opt.wr_max    = 1.01;
opt.w_ndiv    = 100;

%opt.wr_min    = 1.0;
%opt.wr_max    = 1.3;
%opt.w_ndiv    = 200;

opt = [];
opt.kmax      = 0;
opt.realbasis = 0;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
opt.kmax      = 4;
opt.realbasis = 0;
[H1,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
opt.kmax      = 4;
%opt.realbasis = 1;
opt.skewproj = 1;
[H2,freq]   = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

clf;
figure(1);
plot(freq/1e6/2/pi,log10(abs(H1-H)./abs(H)),'r*-');
hold on;
plot(freq/1e6/2/pi,log10(abs(H2-H)./abs(H)),'b*-');
xlabel('Frequency [MHz]');
ylabel('Relative error');
