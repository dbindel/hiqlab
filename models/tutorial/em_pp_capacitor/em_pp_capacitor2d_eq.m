clear;
qcloseall;
% -- Parameters

param.reps = 160;
param.dense = 1e-6;
param.order = 3;

w0  = 1.917998e+08*2*pi;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('em_pp_capacitor2d.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt)

% -- Extract ids
eno_d = Lua_get_double(L,'eno_d')+1;
eno_s = Lua_get_double(L,'eno_s')+1;
idg_d = Lua_get_double(L,'idg_d')+1;
idg_s = Lua_get_double(L,'idg_s')+1;
param.eno = [eno_d,eno_s];
param.idg = [idg_d,idg_s];

% -- Extract LRCC parameters
[eL,eR,eC,eC0,w] = equiv_LRCC(mesh, w0, param);

fprintf('Equivalent circuit parameters:\n');
fprintf('L  = %g\n', eL);
fprintf('R  = %g\n', eR);
fprintf('C  = %g\n', eC);
fprintf('C0 = %g\n', eC0); 

% -- Clean up
Mesh_delete(mesh);
