% HiQLab
% Copyright (c): Regents of the University of California
% $Id: diagnose_pml.m,v 1.1 2006/05/01 04:47:52 tkoyama Exp $

function [rcoeff, kcoeff, resid] = diagnose_pml(param, k)

% -- Compute PML solution

mesh  = Mesh_load('test1d.lua', param);
Mesh_make_harmonic(mesh, k);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
u     = -(K-k^2*M)\F;
Mesh_delete(mesh);

% -- Compute incoming and outgoing solutions

order = param.order;
H   = K-k^2*M;
H11 = full( H(1:order,1:order)         );
H12 = full( H(1:order,order+(1:order)) );
II  = eye(order);
ZZ  = zeros(order);

[V,D] = eig([-H12', -H11; ZZ, II], [ZZ, H12; II, ZZ]);
Vall = V;
Dall = D;
d = diag(D);

[d,I] = sort(d);
V = V(:,I((order-1)+(1:2)));
d = d(I((order-1)+(1:2)));

% -- Decompose PML solution

uu = u(order+(1:2*order));
xx = V\uu;

% -- Print or return diagnostics

rcoeff = abs(xx(2)/xx(1));
kcoeff = k - abs(log(d(1)));
resid  = norm(V*xx-uu);

if nargout == 0
  fprintf('Incoming / outgoing:        %e\n', rcoeff);
  fprintf('Error in wave num:          %e\n', kcoeff);
  fprintf('Discrete solution mismatch: %e\n', resid);
end
