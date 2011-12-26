clear
qcloseall;
% -- PML parameter
param.f0    =20;
param.order = 3;
param.dense = 1e-6;

% -- Load mesh from Lua input file
mesh  = Mesh_load('diskmesh.lua',param);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
w         = 2.889953e+08*2*pi;
drive_pat = Mesh_get_vector(mesh,'bode_drive_func2');
harmonic_state(mesh,drive_pat,w);

% -- Plot forced response
figure(1);
opt.deform  = 1e-6;
opt.axequal = 1;
opt.cfields = [1,2];
plotcycle2d(mesh,opt.deform,opt);
