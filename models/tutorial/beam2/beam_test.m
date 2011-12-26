function beam_test(which_file)

qcloseall;

% -- Load mesh from Lua input file
if nargin < 1, which_file = 1; end
if which_file == 1,  fname = 'beam_test.lua';            end
if which_file == 2,  fname = 'beam_test_blocks2dn.lua';  end
if which_file == 3,  fname = 'beam_test_blocks2d.lua';   end
mesh = Mesh_load(fname);

% -- Display the mesh
qfigure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Find the static displacement: U = K\F
static_state(mesh);

% -- Display displacements
qfigure(2);
opt.deform = 5;
plotfield2d(mesh,opt);

% -- Clean up
Mesh_delete(mesh);
