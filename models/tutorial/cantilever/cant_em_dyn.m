function cant_em_dyn()

params.order = 1;
params.dense = 1e-7;
params.V     = 100;
[mesh, L] = Mesh_load('cant_em.lua', params);

state_opt.nonlinear = 'NR';
static_state(mesh, state_opt);

w0      = 31.018e6*2*pi;
[V,w,Q] = emmode(mesh,w0,1);

qcloseall;
popt.deform  = 1e-6; 
popt.cfields = 3;
popt.cbias   = 1;
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plot_mode(mesh,V,w,popt);

Mesh_delete(mesh);
