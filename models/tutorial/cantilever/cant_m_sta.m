function cant_m_sta(which_file)

if nargin < 1, which_file = 1; end
if which_file == 1,  filename = 'cant_m.lua';          end
if which_file == 2,  filename = 'cant_m_nondim.lua';   end

% -- Load and analyze
[mesh, L] = Mesh_load(filename);
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
