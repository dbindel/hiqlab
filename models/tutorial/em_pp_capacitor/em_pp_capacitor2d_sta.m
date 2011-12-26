clear;
qcloseall;
% -- Parameters
param.reps = 160;
param.dense= 1e-6;
param.order= 3;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('em_pp_capacitor2d.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute for static state
sopt.nonlinear = 'NR';
static_state(mesh,sopt)

% -- Plot voltage solution
figure(2)
popt.cbias   = 1;
popt.cfields = 3;
popt.axequal = 1;
plotfield2d(mesh,popt);

% -- Analytical solution
eps     = Mesh_get_scale(mesh,'eps');
reps    = Lua_get_double(L,'reps');
Cw      = Lua_get_double(L,'Cw');
Ch      = Lua_get_double(L,'Ch');
Ct      = Lua_get_double(L,'Ct');
Can     = reps*eps*Cw*Ct/Ch;

% -- Extract values if electrodes are present
U       = Mesh_get_u(mesh);
Ct      = Lua_get_double(L,'Ct');
ele_t_V = Mesh_get_vector(mesh,'electrode_top_V2');
ele_t_Q = Mesh_get_vector(mesh,'electrode_top_Q2');
Vcomp   = ele_t_V'*U;
Ccomp   = ele_t_Q'*U/Vcomp*Ct;

% -- Print results
fprintf('Analytic capacitance:%d\n',Can);
fprintf('Computed capacitance:%d\n',Ccomp);

% -- Clean up
Mesh_delete(mesh);
