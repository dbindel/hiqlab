% [V,w,Q] = emcmode(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
%
%
% opt must contain following parameters
%   eno - array of element numbers for electrodes
%   idg - array of global number for driving electrode
%
% opt   contains optional parameters:
%   mech - Conduct purely mechanical analysis
%   use_matlab - use matlab eigs or not (Default:1)

function [V,w,Q] = emcmode(mesh, w0, nev, opt);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

mech           = qoptdefault(opt, 'mech',      0);
opt.use_matlab = qoptdefault(opt, 'use_matlab', 1);

if mech
   [V,w,Q] = emcmode_mech_1(mesh,w0,nev,opt);
   return
end

if opt.use_matlab
  [V,w,Q] = emcmode_matlab_1(mesh, w0, nev, opt);
end

% -------------------------------------------------------------------
% [V,w,Q] = emcmode_matlab_1(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values. Right eigenvectors are computed
%
% param must contain following parameters
%   eno - array of element numbers for electrodes
%   idg - array of global number for driving electrode
%
function [V,w,Q] = emcmode_matlab_1(mesh, w0, nev, param);

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

idg            = param.idg;
eno            = param.eno;

numid   = Mesh_get_numid(mesh);    % Number of IDs

% -- Extract nodal ids
numid  = Mesh_get_numid(mesh);
id     = Mesh_get_id(mesh);
ndm    = Mesh_get_ndm(mesh);
ndf    = Mesh_get_ndf(mesh);
id_m   = [];
for i = 1:ndf-2
    tmp = find(id(i,:)>0);
    id_m= union(id_m,id(i,tmp));
end
id_p  = find(id(ndf-1,:)>0);
id_p  = id(ndf-1,id_p);
id_v  = find(id(ndf,:)>0);
id_v  = id(ndf,id_v);

% -- Extract branch and global variable ids of electrode
for i = 1:length(eno)
    id_bq(i) = Mesh_branchid(mesh,1,eno(i));
    id_gv(i) = Mesh_globalid(mesh,idg(i));
end

% -- SHOULD incorporate case when there are other global dofs

% -- Sort ids
id_1 = [id_m,id_p];
id_2 = setdiff([1:numid],[id_1,id_gv,id_bq]);
id_3 = [id_gv,id_bq];
id_12= [id_1,id_2];

nid_1 = length(id_1);
nid_2 = length(id_2);
nid_3 = length(id_3);

% -- Extract matrices
[M,K,C, su,sf] = Mesh_assemble_mkc_nd(mesh);

% -- Nondimensionalization
cT = qoptdefault(param, 'cT', Mesh_get_scale(mesh, 'T'));
M = M / cT^2;
C = C / cT;
w0 = w0*cT;

At      = -K(id_3,id_3)\K(id_3,id_12);
M       = M(id_12,id_12);
C       = C(id_12,id_12) + C(id_12,id_3)*At;
K       = K(id_12,id_12) + K(id_12,id_3)*At;

I11     = speye(nid_1);
Z11     = spalloc(nid_1,nid_1 ,0);
Z12     = spalloc(nid_1,nid_2 ,0);
Z22     = spalloc(nid_2,nid_2 ,0);

% -- Extract system submatrices
I1  = 1:nid_1;
I2  = nid_1+1:(nid_1+nid_2);
K11 = K(I1,I1);  % Device stiffness
K12 = K(I1,I2);  % Electromechanical coupling
K22 = K(I2,I2);  % Conductance matrix
M11 = M(I1,I1);  % Device mass
C21 = C(I2,I1);  % Mechanical to electro coupling
C22 = C(I2,I2);  % Capacitance matrix

% -- Symmetrize matrices that should be symmetric
K11 = (K11 + K11.')/2;
M11 = (M11 + M11.')/2;
C22 = (C22 + C22.')/2;
K22 = (K22 + K22.')/2;

% -- Form generalized eigenvalue problem and solve
B    = [  I11 ,  Z11 ,  Z12 ;
          Z11 ,  M11 ,  Z12 ;
          Z12',  Z12',  C22 ];

A    = [  Z11 ,  I11 ,  Z12 ;
         -K11 ,  Z11 , -K12 ;
          Z12', -C21 , -K22 ];

pmlopt.linearized = 1;
pmlopt.use_matlab = 1;
[VV,w,Q] = pml_mode(B,A,1i*w0,nev,pmlopt);
w = w/1i;
Q = abs(w)./(2*imag(w));
V = zeros(numid,nev);

% -- Shape mode vectors
V(id_12,:) = VV([1:nid_1,2*nid_1+1:2*nid_1+nid_2],:);
V(id_3, :) = At*V(id_12,:);

% -- Redimensionalization
V = spdiag(su) * V;
w = w / cT;


% -------------------------------------------------------------------
% [V,w,Q] = emcmode_mech_1(mesh, w0, nev, opt);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (thermal part is zero) Compute only the mechanical frequency

function [V,w,Q] = emcmode_mech_1(mesh, w0, nev, opt);

if nargin < 3, nev = 1;    end
if nargin < 4, opt = []; end

numid    = Mesh_get_numid(mesh);      % Number of IDs
id       = Mesh_get_id(mesh);
ndf      = Mesh_get_ndf(mesh);
id_m     = [];
for i = 1:ndf-2
    tmp = find(id(i,:)>0);
    id_m= union(id_m,id(i,tmp));
end

nm            = length(id_m);
[Muu,Kuu,su]  = Mesh_assemble_mk_nd(mesh);
Muu           = Muu(id_m,id_m);
Kuu           = Kuu(id_m,id_m);

[VV,w,Q]  = pml_mode(Muu,Kuu,w0,nev);
V         = zeros(numid,nev);
V(id_m,:) = VV;
V = spdiag(su) * V;
