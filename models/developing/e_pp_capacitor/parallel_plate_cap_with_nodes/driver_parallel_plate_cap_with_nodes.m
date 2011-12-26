% Parallel plate capacitance test calculation
% Using global shape functions
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_parallel_plate_cap_with_nodes.m,v 1.1 2006/05/01 07:56:36 tkoyama Exp $

%[mesh,L] = Mesh_load('parallel_plate_cap_with_nodes.lua');
[mesh,L] = Mesh_load('parallel_plate_cap_with_nodes1.lua');
numid    = Mesh_get_numid(mesh);
idg1     = Lua_get_double(L,'idg1')+1;   % Get global variable number of Electrode
eno1     = Lua_get_double(L,'eno1')+1;   % Get element number of Electrode
vid      = Mesh_globalid(mesh,idg1);     % Get ID of global variable
qid      = Mesh_branchid(mesh,1,eno1);   % Get ID of lagrange multiplier
fprintf('Numid:%d\n',numid);

% Plot mesh of capacitor
figure(1);
plotmesh(mesh);
axis equal;

% Find the corresponding voltage distribution
K = Mesh_assemble_dR(mesh, 1,0,0);
F =-Mesh_assemble_R (mesh);
U = K\F;
Mesh_set_u(mesh,U);

% Plot voltage distribution inside capacitor
figure(2);
plotopt.cfields = 1;
plotopt.ufields = 1;
plotfield2d(mesh,plotopt);
axis equal;

% Compute the capacitance as Q/V and by the hand formula
eps0 = Lua_get_double(L, 'eps0');
w    = Lua_get_double(L, 'w'   );
g    = Lua_get_double(L, 'g'   );
C1   = eps0*w/g;
C2   = U(qid)/U(vid);

fprintf('Capacitance by hand:   %e\n', C1);
fprintf('Capacitance from mesh: %e\n', C2);

Mesh_delete(mesh);
