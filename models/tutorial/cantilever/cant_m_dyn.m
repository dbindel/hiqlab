function cant_m_dyn(which_file)

if nargin < 1, which_file = 1; end
if which_file == 1,  filename = 'cant_m.lua';         end
if which_file == 2,  filename = 'cant_m_nondim.lua';  end

[mesh, L] = Mesh_load(filename);

w0      = 0;
[V,w,Q] = mechmode(mesh,w0,1);

qcloseall;
popt.deform  = 1e-6; 
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plot_mode(mesh,V,w,popt);

Mesh_delete(mesh);
