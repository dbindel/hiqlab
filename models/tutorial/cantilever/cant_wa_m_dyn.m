function cant_wa_m_dyn()

params.f0 = 20;
[mesh, L] = Mesh_load('cant_wa_m.lua', params);

w0      = 31.018e6*2*pi;
[V,w,Q] = mechmode(mesh,w0,1);

qcloseall;
popt.deform  = 1e-6; 
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plot_mode(mesh,V,w,popt);

Mesh_delete(mesh);
