clear;
qcloseall;
% -- Pass parameters into Lua file
params.num_elem = 1;
params.order    = 1;

% -- Load mesh from Lua input file
mesh = Mesh_load('circular_disk.lua',params);

% -- Display the mesh
qfigure(1);
opt.axequal = 1;
plotmesh(mesh,opt);
Mesh_get_numid(mesh)

% -- Construct stiffness matrix and forcing vector
%    and solve for displacements
%    U = K\F;
%static_state(mesh);

% -- Display displacements
%figure(2);
%opt.deform = 1;
%plotfield2d(mesh,opt);

% -- Clean up
Mesh_delete(mesh);
