% [wall, Qall] = driver_sweep(w0, nmodes, pname, prange, param, opt)
%
% Return an array of complex frequencies (wall) and Q values (Qall)
% for the nmodes modes closest to w0 in a disk resonator as one
% parameter is varied.  The param structure holds the nominal
% parameters; the one named pname is varied over prange.  opt is
% used to pass options to plotcycle2d if a movie is desired.

function [wall, Qall] = driver_sweep2(w0, nmodes, pname, prange, param, opt)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_sweep_pert.m,v 1.1 2006/05/01 05:06:55 tkoyama Exp $

if nargin < 5,  param = [];  end
if nargin < 6,  opt = [];    end

wall = zeros(nmodes, length(prange));
Qall = zeros(nmodes, length(prange));

if ~isfield(opt, 'sframe')
  opt.sframe = 1; 
end

VV = [];
t0 = clock;
for kk = 1:length(prange);

  param = setfield(param, pname, prange(kk));
  mesh = Mesh_load('diskmesh.lua', param);
  [M,K] = Mesh_assemble_mk(mesh);
  Mesh_delete(mesh);

  nmode = 2;
  if size(VV,1) ~= size(M,1)
    [VV,w,Q] = pml_mode(M,K,w0,nmode, opt);
  else
    Mk = VV.'*M*VV;
    Kk = VV.'*K*VV;
    w = sqrt(eig(Mk\Kk));
    [tmp,I] = sort(abs(w-w0));
    w = w(I);
    Q = abs(w)./(2*imag(w));
  end
  
  % -- Print diagnostic information
  t1 = clock;
  telapsed = etime(t1,t0);
  fprintf('Thickness: %g\n', prange(kk));
  for jj = 1:nmodes
    fprintf('Mode     : %g MHz\t%g Q\n', real(w(jj))/2e6/pi, Q(jj));
  end
  fprintf('Time elapsed   : %g\n', telapsed);
  fprintf('Time remaining : %g\n', telapsed*(length(prange)/kk-1));

  % -- Save diagnostic information
  wall(:,kk) = w;
  Qall(:,kk) = Q;
  
  [Qmax,I] = max(Q);
  wmax = w(I);

  % -- Plot if desired
  if isfield(opt, 'framepng')
    driver_force(real(wmax), param, opt);
    opt.sframe = opt.sframe + opt.nframes;
  end

end
