clear;
qcloseall;
% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('ladder_filter.lua');

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Solve AC case(Transfer function)
pL1 = Lua_get_double(L,'L1');
pC1 = Lua_get_double(L,'C1');
Rs  = Lua_get_double(L,'Rs');
Rl  = Lua_get_double(L,'Rl');
wc  = 1/sqrt(pC1*pL1);

% -- Get forcing and sensing patterns
sense_patA = Mesh_get_vector(mesh, 'bode_senseA');
sense_patB = Mesh_get_vector(mesh, 'bode_senseB');
drive_pat  = Mesh_get_vector(mesh, 'ac_source');
sense_pat  = [sense_patA,sense_patB];

opt.wr_min = 0.95;
opt.wr_max = 1.05;
opt.w_ndiv = 100;
opt.mkc    =  1;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
H = H(2,:)./H(1,:)*(Rs+Rl)/Rl;

% -- Compute analytical result
[Ha] = ladder_filter_analytical(L,freq);

% -- Show bode plot
figure(2);
opt.lstyle    = 'r';
plot_bode(freq,H,opt);
hold on;
opt.lstyle    = 'bs';
plot_bode(freq,Ha,opt);
hold off;

% -- Clean up
Mesh_delete(mesh);
