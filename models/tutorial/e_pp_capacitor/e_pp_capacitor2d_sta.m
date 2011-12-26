function e_pp_capacitor2d_sta(which_version)

qcloseall;

% -- Parameters
param.reps = 160;
param.dense= 1e-6;

if nargin < 1, which_version = 1; end
if which_version == 0, param.use_case = ''; end
if which_version == 1, param.use_case = 'eppc2d_using_globals_sta.lua'; end
if which_version == 2, param.use_case = 'eppc2d_using_electrodes_sta.lua'; end

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('e_pp_capacitor2d.lua',param);

% -- Display the mesh
qfigure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute for static state (f0 = 0)
static_state(mesh)

% -- Plot voltage solution
qfigure(2)
popt.cbias   = 1;
popt.ufields = 1;
popt.cfields = 1;
popt.axequal = 1;
plotfield2d(mesh,popt);

% -- Analytical solution
eps     = Mesh_get_scale(mesh,'eps');
reps    = Lua_get_double(L,'reps');
Cw      = Lua_get_double(L,'Cw');
Ch      = Lua_get_double(L,'Ch');
Ct      = Lua_get_double(L,'Ct');
Can     = reps*eps*Cw*Ct/Ch;

% -- Find top and bottom charge, compute solution
force1  = Mesh_get_force(mesh);
top_Q   = Mesh_get_vector(mesh,'sense_top_Q2',2);
bot_Q   = Mesh_get_vector(mesh,'sense_bot_Q2',2);
Qtop    = sum(sum(top_Q.*force1));
Qbot    = sum(sum(bot_Q.*force1));
Vdc     = Lua_get_double(L,'Vf_dc');
cL      = Mesh_get_scale(mesh,'L');
Ct      = Lua_get_double(L,'Ct');
Ccomp   = Qtop/Vdc*(Ct/cL);

% -- Print results
fprintf('Analytic capacitance:%e\n',Can);
fprintf('Computed capacitance:%e\n',Ccomp);

% -- Extract values if electrodes are present
if strcmp(param.use_case,'eppc2d_using_electrodes_sta.lua')
    U = Mesh_get_u(mesh);

    % -- Method 1
    ele_t_V = Mesh_get_vector(mesh,'electrode_top_V2');
    ele_t_Q = Mesh_get_vector(mesh,'electrode_top_Q2');
    Vm1 = ele_t_V'*U;
    Cm1 = ele_t_Q'*U/Vm1*(Ct/cL);

    % -- Method 2
    ele_b_V = Mesh_get_vector(mesh,'electrode_bot_V2');
    ele_b_Q = Mesh_get_vector(mesh,'electrode_bot_Q2');
    Vm2 = ele_b_V'*U;
    Cm2 =-ele_b_Q'*U/Vm1*(Ct/cL);

    fprintf('Computed cap(glob_t):%e\n',Cm1);
    fprintf('Computed cap(glob_b):%e\n',Cm2);
end

% -- Clean up
Mesh_delete(mesh);
