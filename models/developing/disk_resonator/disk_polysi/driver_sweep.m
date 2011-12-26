% [wall, Qall] = driver_sweep(w0, nmodes, pname, prange, param, opt)
%
% Return an array of complex frequencies (wall) and Q values (Qall)
% for the nmodes modes closest to w0 in a disk resonator as one
% parameter is varied.  The param structure holds the nominal
% parameters; the one named pname is varied over prange.  opt is
% used to pass options to plotcycle2d if a movie is desired.

function [wall, Qall] = driver_sweep(w0, nmodes, pname, prange, param, opt)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_sweep.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $

if nargin < 5,  param = [];  end
if nargin < 6,  opt = [];    end

wall = zeros(nmodes, length(prange));
Qall = zeros(nmodes, length(prange));

if ~isfield(opt, 'sframe')
  opt.sframe = 1; 
end

t0 = clock;
for kk = 1:length(prange);
 
  param = setfield(param, pname, prange(kk));
  [w,Q] = driver_mode(w0, param, nmodes);

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
