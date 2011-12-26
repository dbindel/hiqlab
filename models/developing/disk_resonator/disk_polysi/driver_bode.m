% [H] = driver_bode(freq, params, kmax, w0)
%
% Compute a Bode plot using the given frequency sampling points.
%
% If kmax and w0 are given, use model reduction via an Arnoldi
% expansion about the shift w0.  If the shift frequency w0 is omitted,
% use w0 = mean(freq).

function [H] = driver_bode(freq, params, kmax, w0)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_bode.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $

if nargin < 2, params = [];     end
if nargin < 3, kmax = 0;        end
if nargin < 4, w0 = mean(freq); end

mesh  = Mesh_load('diskmesh.lua', params);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
L     = F / sum(F);
Mesh_delete(mesh);

if kmax > 0
  opt.realbasis = 1;
  [M,K,L,F] = rom_arnoldi(M,K,L,F, kmax, w0, opt);
end

use_umfpack = ((kmax == 0) & exist('umfpack'));
H = zeros(length(freq),1);

for idx = 1:length(freq)

  w = freq(idx);
  if (kmax == 0)
    fprintf('%d: %d Hz\n', idx, freq(idx)); 
  end

  if use_umfpack,   H(idx) = L' * umfpack(K - w^2*M, '\', F); 
  else              H(idx) = L' * ((K - w^2*M)        \   F);
  end

end

figure(1);
plot_bode(freq, H);

