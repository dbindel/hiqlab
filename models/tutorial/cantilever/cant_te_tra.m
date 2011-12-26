function cant_te_tra()

[mesh, L] = Mesh_load('cant_te.lua');

wc          = 31.018e6*2*pi;
bopt.wr_min = 0.99;
bopt.wr_max = 1.01;
bopt.w_ndiv = 100;
bopt.kmax   = 0;
bopt.mkc    = 1;

drive_pat = Mesh_get_vector(mesh,'bode_force_function2');
sense_pat = Mesh_get_vector(mesh,'bode_sense_function2');
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,bopt);

qcloseall;
opt.axequal = 1;
qfigure(1); plotmesh(mesh,opt);
figure(2);  plot_bode(freq,H);

Mesh_delete(mesh);
