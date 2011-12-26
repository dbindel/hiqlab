% [V,w,Q] = mechmode(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
%
% param contains optional parameters:
%   use_matlab - use matlab eigs or not (Default:0)

function [V,w,Q] = mechmode(mesh, w0, nev, param);

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

opt.use_matlab = qoptdefault(param, 'use_matlab', 0);

[M,K, su]  = Mesh_assemble_mk_nd(mesh);
[V,w,Q] = pml_mode(M,K,w0,nev,opt);
V = spdiag(su)*V;

V = V / diag(sqrt(diag(V'*V)));
Q = abs(w)./(2*imag(w));
