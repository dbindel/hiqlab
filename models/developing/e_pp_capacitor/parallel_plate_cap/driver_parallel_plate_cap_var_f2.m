% Parallel plate capacitance test calculation(variant)
% Uses TieField2 Element to tie
%  1. Bottom potentials
%  2. Top left potentials
%  3. Top right potentials
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_parallel_plate_cap_var_f2.m,v 1.1 2006/05/01 07:55:57 tkoyama Exp $

[mesh,L] = Mesh_load('parallel_plate_cap_var_f2.lua');
numid = Mesh_get_numid(mesh);
fprintf('Numid:%d\n',numid);

% Plot mesh of capacitor
figure(1);
plotmesh(mesh);
axis equal;

% Obtain info on TieField elements
ntie1  = Lua_get_double(L, 'ntie1')+1;    % Element No.
ntie1I = Mesh_branchid(mesh,1,ntie1);     % ID No. for tied Voltage
ntie1A = Mesh_nbranch_idja(mesh,ntie1);   % No. of active branch variables
                                          % (Since some nodes may be tied)
ntie1E = Mesh_branchid(mesh,ntie1A,ntie1);% ID No. for total charge
                                          % (Last branch variable)
ntie2  = Lua_get_double(L, 'ntie2')+1;
ntie2I = Mesh_branchid(mesh,1,ntie2);
ntie2A = Mesh_nbranch_idja(mesh,ntie2);
ntie2E = Mesh_branchid(mesh,ntie2A,ntie2);
ntie3  = Lua_get_double(L, 'ntie3')+1;
ntie3I = Mesh_branchid(mesh,1,ntie3);
ntie3A = Mesh_nbranch_idja(mesh,ntie3);
ntie3E = Mesh_branchid(mesh,ntie3A,ntie3);

% Set 1. Bottom    potential to -2
%     2. Top left  potential to  1
%     3. Top right potential to  3
qq    = zeros(numid,1);
qq(ntie1E) = -2;
qq(ntie2E) =  1;
qq(ntie3E) =  3;

% Find the corresponding voltage distribution
K  = Mesh_assemble_dR(mesh, 1,0,0);
vv = K\qq;

% Plot voltage distribution inside capacitor
Mesh_set_u(mesh,vv);
plotopt.ufields = 1;
plotopt.cfields = 1;
plotfield2d(mesh,plotopt);
axis equal;

Mesh_delete(mesh);
