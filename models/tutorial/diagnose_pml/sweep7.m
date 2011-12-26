% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep7.m,v 1.1 2006/05/01 04:47:50 tkoyama Exp $

% Vary damping per wavelength and elements per wavelength

param.epw   = 10;
param.dpw   = 40;
param.Ne1   = 10;
param.Ne2   = 150;
param.order = 1;
k           = 2*pi/param.epw;

rl = [];
epwl = 3:0.5:30;
dpel = logspace(-2,1,50);

for idpe = 1:length(dpel)
  for iepw = 1:length(epwl)
    param.dpw = dpel(idpe) * epwl(iepw);
    param.epw = epwl(iepw);
    k         = 2*pi/param.epw;
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpe,iepw) = rcoeff;
    kl(iepw) = kcoeff;
  end
end

figure(1);
[cs,h] = contour(epwl, log10(dpel), log10(rl)); 
xlabel('Elements per wave (kh)^{-1}');
ylabel('Log_{10} (\beta h)');
clabel(cs,h);
