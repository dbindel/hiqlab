clear;
qcloseall;
% -- Parameters
param.use_case  = 'using_globals_sta.lua';
%param.use_case  = 'using_electrodes_sta.lua';

param.order        = 1;
param.dense        = 1.0e-6;

param.fill_gap     = 1;
param.reps         = 160;
param.Vf_dc        = 100;

param.f0           =   0;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('dielectric_drive.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt);

% -- Check for pull-in
U              = Mesh_get_u(mesh);
sense_gap_dist = Mesh_get_sense_u(mesh,'sense_gap_dist');
Dt             = Lua_get_double(L,'Dt');
fprintf('Gap distance / 3     [um]:%d\n',               Dt/3/1e-6);
fprintf('Gap distance decrease[um]:%d\n',(sense_gap_dist'*U)/1e-6);

% -- Plot voltage solution
figure(2)
popt.cbias   = 1;
popt.cfields = 3;
popt.axequal = 1;
plotfield2d(mesh,popt);

% -- Extract values
% -- Find top and bottom charge, compute solution
force1  = Mesh_get_force(mesh);
top_Q   = Mesh_get_sense_force(mesh,'sense_top_Q');
bot_Q   = Mesh_get_sense_force(mesh,'sense_bot_Q');
Qtop    = sum(sum(top_Q.*force1));
Qbot    = sum(sum(bot_Q.*force1));
Vdc     = Lua_get_double(L,'Vf_dc');
cL      = Mesh_get_scale(mesh,'L');
Ct      = Lua_get_double(L,'Lt');
Ccomp   = Qtop*(Ct/cL)/Vdc;

% -- Analytical solution
eps     = Mesh_get_scale(mesh,'eps');
reps    = Lua_get_double(L,'reps');
Cw      = Lua_get_double(L,'Bw')*2;
Ch      = Lua_get_double(L,'Dt');
Ct      = Lua_get_double(L,'Lt');
Can     = reps*eps*Cw*Ct/Ch;

% -- Print results
fprintf('Computed capacitance:%d\n',Ccomp);
fprintf('Analytic capacitance:%d\n',Can);
