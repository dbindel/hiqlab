% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep3.m,v 1.1 2006/05/01 04:47:54 tkoyama Exp $

% Vary damping per wavelength and elements per wavelength

param.epw   = 10;
param.dpw   = 40;
param.Ne1   = 10;
param.Ne2   = 180;
param.order = 3;
k           = 2*pi/param.epw;

rl = [];
%dpwl = 1:20;
dpwl = logspace(-1,1,20);
epwl = 5:0.5:30;

for idpw = 1:length(dpwl)
  for iepw = 1:length(epwl)
    param.dpw = dpwl(idpw);
    param.epw = epwl(iepw);
    k         = 2*pi/param.epw;
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpw,iepw) = rcoeff;
  end
end

[cs,h] = contour(epwl, log10(dpwl), log10(rl)); 
xlabel('Elements per wavelength');
ylabel('Log_{10} (\beta)');
clabel(cs,h);
