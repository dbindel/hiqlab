function cant_wa_m_sta(which_file)

% -- Load and analyze
params.f0 = 0;
[mesh, L] = Mesh_load('cant_wa_m.lua', params);
static_state(mesh);

% -- Display mesh, displacement fields, and tip displacement
qcloseall;
popt.deform  = 100;
popt.axequal = 1;
qfigure(1); plotmesh(mesh,popt);
qfigure(2); plotfield2d(mesh,popt);
ytip = -Lua_get_double(L,'w')/2;
print_tip(mesh, L, ytip);

Mesh_delete(mesh);
