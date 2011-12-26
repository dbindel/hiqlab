function cant_wa_te_sta()

% -- Load and analyze
params.f0 = 0;
[mesh, L] = Mesh_load('cant_wa_te_sta.lua', params);
static_state(mesh);

% -- Display mesh, displacement fields, and tip displacement
qcloseall;
popt.deform  = 1000;
popt.cfields = 3;
popt.cbias   = 1;
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plotfield2d(mesh,popt);
ytip = -Lua_get_double(L,'w')/2;
print_tip(mesh, L, ytip);

Mesh_delete(mesh);
