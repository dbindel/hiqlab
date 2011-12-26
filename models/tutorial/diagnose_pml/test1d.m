% Driver for simple PML test: block sitting on a half space.
% Version 1: Use Lua for everything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test1d.m,v 1.1 2006/05/01 04:47:58 tkoyama Exp $

% -- Assemble and solve linear system

param.epw   = 20;
param.dpw   = 7;
param.Ne1   = param.epw;
param.Ne2   = 10;
param.order = 3;
k           = 2*pi/param.epw;

mesh  = Mesh_load('test1d.lua', param);
Mesh_make_harmonic(mesh, k);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
u     = -(K-k^2*M)\F;
Mesh_set_u(mesh, u);

% -- Plot forced response

x = Mesh_get_x(mesh);
u = Mesh_get_disp(mesh);
%v = exp(-1i*k*(x+param.Ne1));
v = exp(-1i*k*(x+param.Ne1) - k*(param.dpw/param.epw)*max(0,x).^2/2);

plot(x, real(u), 'b.',   x,  real(v), 'b-', ...
     x, imag(u), 'r.',   x,  imag(v), 'r-');

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

[junk,I] = sort(d-1);
V = V(:,I(1:2));
d = d(I(1:2));

uu = u((order+1)+(1:2*order)).';
xx = V\uu;

fprintf('Discrete solution mismatch: %e\n', norm(V*xx-uu));
fprintf('Fraction error in wave num: %e\n', abs( log(d(1)) + k*1i ));
fprintf('Incoming / outgoing:        %e\n', abs(xx(2)/xx(1)));

