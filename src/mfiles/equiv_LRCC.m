% [eL,eR,eC,eC0,w] = equiv_LRCC(mesh, w0, param);
%
% Compute equivalent LRCC parameters near target frequency w0
%
% param must contain following parameters
%   eno - array of element numbers for electrodes
%   idg - array of global number for driving electrode

function [eL,eR,eC,eC0,w] = equiv_LRCC(mesh, w0, param);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, param = []; end

cT      = qoptdefault(param, 'cT', Mesh_get_scale(mesh, 'T'));
idg     = param.idg;
eno     = param.eno;
numid   = Mesh_get_numid(mesh);    % Number of IDs

w0 = w0 * cT;

% -- Extract nodal ids
numid  = Mesh_get_numid(mesh);
id     = Mesh_get_id(mesh);
ndm    = Mesh_get_ndm(mesh);
ndf    = Mesh_get_ndf(mesh);
id_n   = [];
for i = 1:ndf-1
    tmp = find(id(i,:)>0);
    id_n= union(id_n,id(i,tmp));
end
id_v  = find(id(ndf,:)>0);
id_v  = id(ndf,id_v);

% -- Extract branch and global variable ids
for i = 1:length(eno)
    id_bq(i) = Mesh_branchid(mesh,1,eno(i));
    id_gv(i) = Mesh_globalid(mesh,idg(i));
end

% -- Sort ids
id_d1=[id_n];
id_d2=[id_gv];
id_c =setdiff([1:numid],[id_n,id_v,id_gv,id_bq]);
id_a =[id_bq,id_v,id_c];

% -- Extract matrices
[M,K,D, su,sf] = Mesh_assemble_mkc_nd(mesh);

% -- Nondimensionalize
M = M / cT^2;
D = D / cT;

% -- Extract system submatrices
K11 = K(id_d1,id_d1);
K12 = K(id_d1,id_d2);
K21 = K(id_d2,id_d1);
K22 = K(id_d2,id_d2);
D11 = D(id_d1,id_d1);
M11 = M(id_d1,id_d1);

% -- Symmetrize matrices that should be symmetric
K11 = (K11 + K11.')/2;
M11 = (M11 + M11.')/2;
D11 = (D11 + D11.')/2;
K22 = (K22 + K22.')/2;

% -- Compute eigenvector
pmlopt.linearized = 1;
nid_d1 = length(id_d1);
I11    = speye(nid_d1,nid_d1);
Z11    = spalloc(nid_d1,nid_d1,0);

% -- Compute right eigenvector
B    = [  I11 ,  Z11;
          Z11 ,  M11;];

A    = [  Z11 ,  I11;
         -K11 , -D11;];

[Vr,w,Q] = pml_mode(B,A,1i*w0,1,pmlopt);
w  = w/1i;
Vr = Vr([1:nid_d1],:);
Q = abs(w)./(2*imag(w));

% -- Compute left eigenvector
B    = [  I11 ,  Z11;
          Z11 ,  M11';];

A    = [  Z11 ,  I11;
         -K11', -D11';];

[Vl,w,Q] = pml_mode(B,A,1i*w0,1,pmlopt);
w  = conj(w/1i);
Vl = Vl([1:nid_d1],:);
Q = abs(w)./(2*imag(w));

% -- Offset capacitance
Lts  = Mesh_get_scale(mesh,'Lt')/Mesh_get_scale(mesh,'L');
K22  = K22 - K21*(K11\K12);
qC0  = K22(1,1) * Lts;

% -- Equivalent mechanical parameters
eM   = Vl'*M11*Vr   * Lts;
eD   = Vl'*D11*Vr   * Lts;
eK   = Vl'*K11*Vr   * Lts;
eta_r= Vl'*K12(:,1) * Lts;
eta_l= K21(1,:)*Vr  * Lts;

% -- Equivalent LRCC parameters
eta2 = eta_l*eta_r;
qL  = eM/eta2;
qR  = eD/eta2;
qC  = eta2/eK;

% -- Convert to real parameters
eC0  =  real(qC0);
eL   = -real(qL);
eR   =  imag(qL)*real(w) + imag(qC)/abs(qC)^2/real(w);
eC   = -abs(qC)^2/real(qC);

% -- Dimensionalize
eL  = eL   * Mesh_get_scale(mesh,'Li');
eR  = eR   * Mesh_get_scale(mesh,'R' );
eC  = eC   * Mesh_get_scale(mesh,'C' );
eC0 = eC0  * Mesh_get_scale(mesh,'C' );
w   = w/cT;

