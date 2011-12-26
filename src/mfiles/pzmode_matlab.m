% [V,w,Q] = pzmode_matlab(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and potential parts).
%
% *Note*: pzmode will reorder the degrees of freedom in the mesh.

function [V,w,Q] = pzmode_matlab(mesh, w0, nev, param);

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

opt             = [];
opt.disp        = qoptdefault(opt, 'disp',        0);
opt.use_umfpack = qoptdefault(opt, 'use_umfpack', exist('umfpack') == 3);
opt.issym       = qoptdefault(opt, 'issym'      , 0);

numid    = Mesh_get_numid(mesh);    % Number of IDs
nm       = pz_block_mesh(mesh);     % Number of mech IDs
np       = numid - nm;              % Number of potential IDs
[M,K,su] = Mesh_assemble_mk_nd(mesh);
M        = (M + M.')/2;
K        = (K + K.')/2;

Im      = 1:nm;
Ip      = (nm+1):(numid);


% -- Check matrices properties
if (isreal(K) & isreal(M)) & isreal(w0)
  opt.isreal = 1;
else
  opt.isreal = 0;
end;

if np == 0  % For the case when there are no electrical field dofs

  [L1,U1,P1,Q1] = qlu(K-w0^2*M);
  [V,D] = eigs(@afunp, length(M), nev, w0^2, opt, L1, U1, P1, Q1, M);
  d = diag(D);

else % np > 0

  [L1,U1,P1,Q1]= qlu(K-w0^2*M);
  V            = zeros(numid,nev);
  [V(Im,:),D]  = eigs(@afun, length(Im), nev, w0^2, opt, L1, U1, P1, Q1, Im, M);
  d            = diag(D);
  V(Ip,:)      =   K(Ip,Im)*V(Im,:);
  V(Ip,:)      = - K(Ip,Ip)\V(Ip,:);
end

V   = spdiag(su)*V;
V   = V  / diag(sqrt(diag(V'*V)));
w   = sqrt(d);
Q   = abs(w)./(2*imag(w));


% ---
% Apply y = (K-sigma*M)\(M*z)
%
function [y] = afun(x, L1, U1, P1, Q1, Im, M)
z     = zeros(length(L1),1);
z(Im) = x;
y1    = Q1*(U1\(L1\(P1*(M*z))));
y     = y1(Im);


% ---
function [y] = afunp(x, L1, U1, P1, Q1, M)
y = Q1*(U1\(L1\(P1*(M*x))));
