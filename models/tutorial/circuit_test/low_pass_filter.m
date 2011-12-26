clear;
qcloseall;
% -- Load mesh from Lua input file
[mesh,L] = Mesh_load('low_pass_filter.lua');

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Solve AC case(Transfer function)
res = Lua_get_double(L,'R');
cap = Lua_get_double(L,'C');
wc  = 1/res/cap;

% -- Get driving and sensing patterns
sense_patA = Mesh_get_vector(mesh, 'bode_senseA');
sense_patB = Mesh_get_vector(mesh, 'bode_senseB');
drive_pat  = Mesh_get_vector(mesh, 'ac_source');
sense_pat = [sense_patA, sense_patB];

% -- Compute transfer function
opt.wr_min = -2;
opt.wr_max =  3;
opt.w_ndiv = 100;
opt.w_type = 'log';
opt.mkc    =  1;
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,opt);
H = H(1,:)./H(2,:);

% -- Show bode plot
figure(2);
opt.logf = 1;
plot_bode(freq/wc,H,opt);

% -- Clean up
Mesh_delete(mesh);
