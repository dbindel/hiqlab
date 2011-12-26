clear;
qcloseall;
% -- Parameters
param.reps     = 160;
param.dense    = 1e-6;
param.order    = 3;

param.mech = 0;
nev        = 1;
w0  = 1.917995e+08*2*pi;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('em_pp_capacitor2d.lua',param);

% -- Compute for static state
sopt.nonlinear = 'NR';
static_state(mesh,sopt)

% -- Extract element & global numbers for electrodes
eno_d = Lua_get_double(L,'eno_d')+1;
eno_s = Lua_get_double(L,'eno_s')+1;
idg_d = Lua_get_double(L,'idg_d')+1;
idg_s = Lua_get_double(L,'idg_s')+1;
param.eno = [eno_d,eno_s];
param.idg = [idg_d,idg_s];

% -- Compute frequencies, Qs, and mode shapes
[V,w,Q] = emcmode(mesh, w0, nev, param);

% -- Show modes
figure(2);
popt.deform  = 1e-5;
popt.axequal = 1;
popt.cfields = [1,2];
popt.cbias   = 10;
popt.nframes = 8;
popt.animate = 1;
plot_mode(mesh,V,w,popt);

% -- Clean up
Mesh_delete(mesh);
