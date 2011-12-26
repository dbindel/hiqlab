function arch

qcloseall;

% -- Pass parameters into Lua file
params.num_elem = 3;
params.order    = 3;

% -- Load mesh from Lua input file
mesh = Mesh_load('arch.lua', params);

% -- Display the mesh
qfigure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Find the static displacement U = K\F
static_state(mesh);

% -- Display displacements
qfigure(2);
opt.deform = 1;
plotfield2d(mesh,opt);

% -- Clean up
Mesh_delete(mesh);
