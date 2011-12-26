% Parallel plate capacitance test calculation
% Using TieField Element
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_parallel_plate_cap_f.m,v 1.1 2006/05/01 07:55:52 tkoyama Exp $

[mesh,L] = Mesh_load('parallel_plate_cap_f.lua');
numid    = Mesh_get_numid(mesh);
fprintf('Numid:%d\n',numid);

% Plot mesh of capacitor
figure(1);
plotmesh(mesh);
axis equal;

% Set the total charge (dual to the tied voltage) to 1
ntie  = Lua_get_double(L, 'ntie')+1; % Get element number of TieField Element
xtie  = Mesh_branchid(mesh,1,ntie);  % Get ID of 1st variable in TieField Element
qq    = zeros(numid,1);
qq(xtie) = 1;

% Find the corresponding voltage distribution
K  = Mesh_assemble_dR(mesh, 1,0,0);
vv = K\qq;
Mesh_set_u(mesh,vv);

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
C2   = qq(xtie)/vv(xtie);

fprintf('Capacitance by hand:   %e\n', C1);
fprintf('Capacitance from mesh: %e\n', C2);

Mesh_delete(mesh);
