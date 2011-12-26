clear;
qcloseall;
% -- Parameters
param.use_case  = 'using_electrodes_sta.lua';

param.order        = 1;
param.dense        = 1.0e-6;

param.fill_gap     = 1;
param.reps         = 160;
param.Vf_dc        = 100;

param.f0           =   0;
f0                 =  20;

w0          = 140.020e6*2*pi;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('dielectric_drive.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt);
U              = Mesh_get_u(mesh);

% -- Reload mesh with PML
Mesh_delete(mesh);
param.f0  =   f0;
[mesh, L] = Mesh_load('dielectric_drive.lua',param);
Mesh_set_u(mesh,U);

% -- Extract ids
eno_d     = Lua_get_double(L,'eno_d')+1;
eno_s     = Lua_get_double(L,'eno_s')+1;
idg_d     = Lua_get_double(L,'idg_d')+1;
idg_s     = Lua_get_double(L,'idg_s')+1;
param.eno = [eno_d,eno_s];
param.idg = [idg_d,idg_s];


% -- Extract LRCC parameters
[eL,eR,eC,eC0,w] = equiv_LRCC(mesh, w0, param);
