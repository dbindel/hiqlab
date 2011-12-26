function cant_em_sta()

% -- Load and analyze
params.order = 1;
params.dense = 1e-7;
params.V     = 100;
[mesh, L] = Mesh_load('cant_em.lua', params);

state_opt.nonlinear = 'NR';
state_opt.R_tol = 1e-14;
static_state(mesh, state_opt);

% -- Display mesh, displacement fields, and tip displacement
qcloseall;
popt.deform  = 10;
popt.cfields = 3;
popt.cbias   = 1;
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plotfield2d(mesh,popt);
ytip = Lua_get_double(L,'g');
print_tip(mesh, L, ytip);

Mesh_delete(mesh);
