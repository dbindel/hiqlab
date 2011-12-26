
% [rcoeff, kcoeff, resid] = diagnose_pml1d(k, order, c, p, N)
%
% Compute the reflection coefficient and dispersion level for a 1D
% problem with a PML of the form
%   lambda(x) = 1 - c*x^p for x > 0
%
% Inputs:
%   k      - desired wave number
%   order  - element order
%   c      - PML profile coefficient
%   p      - PML profile exponent
%   N      - number of elements in the PML
%
% Outputs:
%   rcoeff - discrete reflection coefficient
%   kcoeff - k - kactual, where kactual is the discrete wave number
%   resid  - residual in wave decomposition

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: diagnose_pml1d.m,v 1.1 2005/04/27 17:59:25 dbindel Exp $

function [rcoeff, kcoeff, resid] = diagnose_pml1d(k, order, c, p, Ne2)

% -- Compute PML solution

Ne1  = 10;
mesh = Mesh_new(1);
D    = PMLScalar1d_new(1,1);
Mesh_own_elt(mesh,D);
Mesh_add_block1d(mesh, -Ne1, Ne2, order*(Ne1+Ne2)+1, D, order);
Mesh_initialize(mesh);

pp = Mesh_get_x(mesh);
e  = Mesh_get_e(mesh);
id = Mesh_get_id(mesh);
PMLElement_set_stretch(D, c * ((pp > 0) .* pp).^p);

[M,K]    = Mesh_assemble_mk(mesh);
KK       = K-k^2*M;
R        = KK(2:end-1,1);
KK       = KK(2:end-1,2:end-1);
Mesh_delete(mesh);

ur  = -KK\R;
 
% -- Compute incoming and outgoing solutions

nn    = order;
H11   = full( KK( nn+(1:nn),   nn+(1:nn)) );
H12   = full( KK( nn+(1:nn), 2*nn+(1:nn)) );
II    = eye(nn);
ZZ    = zeros(nn);
  
[V,D] = eig([-H12', -H11; ZZ, II], [ZZ, H12; II, ZZ]);
Vall = V;
Dall = D;
d = diag(D);

[junk,I] = sort(d-1);
V = V(:,I);
d = d(I);
if imag(d(1)) > 0
  d([1;2])   = d([2; 1]);
  V(:,[1,2]) = V(:,[2, 1]);
end

% -- Decompose PML solution

uu = ur(1:2*nn);
xx = V\uu;
%[abs(xx), d]

% -- Print or return diagnostics

rcoeff = abs(xx(2)/xx(1));
kcoeff = k - abs(log(d(1)));
resid  = norm(V*xx-uu);

if nargout == 0
  fprintf('Incoming / outgoing:        %e\n', rcoeff);
  fprintf('Fraction error in wave num: %e\n', kcoeff);
  fprintf('Discrete solution mismatch: %e\n', resid);
end
