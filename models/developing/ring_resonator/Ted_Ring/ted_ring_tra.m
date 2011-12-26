% -- Paramters
param.f0    = 0;
param.order = 3;
param.dense = 5e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('ted_ring.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute frequencies, Qs, and eigenvectors
ted_block_mesh(mesh);
wc        = 617.58e6*2*pi;
drive_pat = Mesh_get_drive_f(mesh,'bode_force_function');
sense_pat = Mesh_get_sense_u(mesh,'bode_sense_function');

opt.wr_min = 0.990;
opt.wr_max = 1.010;
opt.w_ndiv = 600;
opt.kmax   = 30;
opt.mkc    = 1;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);

% -- Show bode plot
figure(2);
plot_bode(freq,H);
