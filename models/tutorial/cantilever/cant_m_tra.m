function cant_m_tra(which_file)

if nargin < 1, which_file = 1; end
if which_file == 1,  filename = 'cant_m.lua';         end
if which_file == 2,  filename = 'cant_m_nondim.lua';  end

[mesh, L] = Mesh_load(filename);

wc          = 31.018e6*2*pi;
bopt.wr_min = 0.99;
bopt.wr_max = 1.01;
bopt.w_ndiv = 100;

drive_pat = Mesh_get_vector(mesh,'bode_force_function2');
sense_pat = Mesh_get_vector(mesh,'bode_sense_function2');
[H,freq] = second_order_bode(mesh,wc,drive_pat,sense_pat,bopt);

qcloseall;
opt.axequal = 1;
qfigure(1); plotmesh(mesh,opt);
figure(2);  plot_bode(freq,H);

Mesh_delete(mesh);
