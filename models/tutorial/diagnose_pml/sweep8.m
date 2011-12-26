% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep8.m,v 1.1 2006/05/01 04:47:57 tkoyama Exp $

% Vary damping per wavelength and PML length

param.epw   = 10;
param.dpw   = 40;
param.Ne1   = 10;
param.Ne2   = 150;
param.order = 2;
k           = 2*pi/param.epw;

rl = [];
nel  = 1:30;
dpel = logspace(-2,1,50);

for idpe = 1:length(dpel)
  param.Ne2 = 150;
  param.dpw = dpel(idpe) * param.epw;
  k         = 2*pi/param.epw;
  [rcoeff_pml, kcoeff, resid] = diagnose_pml(param, k);
  for ine = 1:length(nel)
    param.Ne2 = nel(ine);
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpe,ine)         = rcoeff;
    rl_pml(idpe,ine)     = rcoeff_pml;
    rl_nominal(idpe,ine) = exp(-dpel(idpe) * k * nel(ine)^2);
  end
end

set(0, 'defaulttextinterpreter', 'none');

figure(1);
[cs,h] = contour(nel, log10(dpel), -log10(rl), 1:4); 
xlabel('Number of PML elements');
ylabel('$\log_{10}(\beta h)$');
title('$-\log_{10}(r)$ at $(kh)^{-1} = 10$');
clabel(cs,h);
laprint(1, 'sweep8_1', 'width', 6);

figure(2);
[cs,h] = contour(nel, log10(dpel), -log10(rl_nominal + rl_pml), 1:4); 
xlabel('Number of PML elements');
ylabel('$\log_{10}(\beta h)$');
title('$-\log_{10}(r_{\mathrm{interface}} + r_{\mathrm{nominal}})$ at $(kh)^{-1} = 10$');
clabel(cs,h);
laprint(2, 'sweep8_2', 'width', 6);
