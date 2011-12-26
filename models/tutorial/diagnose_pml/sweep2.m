% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep2.m,v 1.1 2006/05/01 04:47:47 tkoyama Exp $

% Vary damping per element and number of elements

param.epw   = 30;
param.dpw   = 40;
param.Ne1   = 10;
param.Ne2   = 100;
param.order = 3;
k           = 2*pi/param.epw;

rl   = [];
dpel = logspace(-2,2,100);
nel  = 1:30;

for idpe = 1:length(dpel)
  for ine = 1:length(nel)
    param.dpw = dpel(idpe) * param.epw;
    param.Ne2 = nel(ine);
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpe,ine) = rcoeff;
  end
end

figure(1);
[cs,h] = contour(nel, log10(dpel), log10(rl)); 
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
clabel(cs,h);

figure(2);
surf(nel, log10(dpel), log10(rl));
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
zlabel('Log10(Reflection coeff)');
