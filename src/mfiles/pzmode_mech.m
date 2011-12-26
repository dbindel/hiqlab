% [V,w,Q] = pzmode_mech(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and potential parts).
%
% *Note*: pzmode will reorder the degrees of freedom in the mesh.

function [V,w,Q] = pzmode_mech(mesh, w0, nev, param);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

numid   = Mesh_get_numid(mesh);    % Number of IDs
nm      = pz_block_mesh(mesh);     % Number of mech IDs
nt      = numid - nm;              % Number of potential IDs
[M,K] = Mesh_assemble_mk(mesh);

Iu = 1:nm;
Ip = nm+1:numid;

% -- Extract system submatrices
Kuu = K(Iu,Iu);  % Mechanical stiffness
Muu = M(Iu,Iu);  % Mechanical mass

% -- Symmetrize matrices that should be symmetric
Kuu = (Kuu + Kuu.')/2;
Muu = (Muu + Muu.')/2;
[V,w,Q] = pml_mode(Muu,Kuu,w0,nev);
