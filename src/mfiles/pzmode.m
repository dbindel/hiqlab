% [V,w,Q] = pzmode(mesh, w0, nev, param);
%
% Compute nev complex frequencies w (in rad/s) near target frequency
% w0, and also associated Q values.  V contains the eigenvectors
% (mechanical and potential parts).
%
% *Note*: pzmode will reorder the degrees of freedom in the mesh.

function [V,w,Q] = pzmode(mesh, w0, nev, param);

if nargin < 3, nev = 1;    end
if nargin < 4, param = []; end

opt.use_matlab = qoptdefault(param, 'use_matlab', 0);

if opt.use_matlab
  [V,w,Q] = pzmode_matlab(mesh, w0, nev, param);
else
  numid   = Mesh_get_numid(mesh);    % Number of IDs
  nm      = pz_block_mesh(mesh);     % Number of mech IDs
  nt      = numid - nm;              % Number of thermal IDs
  [d,V]   = pz_mode_2(mesh, w0, nev);
  V       = V / diag(sqrt(diag(V'*V)));
  Q       = abs(w)./(2*imag(w));
end
