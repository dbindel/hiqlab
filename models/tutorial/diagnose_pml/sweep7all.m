% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep7all.m,v 1.1 2006/05/01 04:47:56 tkoyama Exp $

% Vary damping per wavelength and elements per wavelength

param.Ne1   = 10;
param.Ne2   = 150;
set(0, 'defaulttextinterpreter', 'none');

eltname = {'linear elements', ...
	   'quadratic elements', ...
	   'cubic elements'};

for p = 1:3
  
  param.order = p;
  
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
  
  figure(p);
  [cs,h] = contour(epwl, log10(dpel), -log10(rl)); 
  clabel(cs,h);
  
  xlabel('Elements per wave $(kh)^{-1}$');
  ylabel('$\log_{10}(\beta h)$');
  title(['$-\log_{10}(r_{\mathrm{interface}})$: ', eltname{p}]);
  laprint(p, sprintf('sweep7_%d', p), 'width', 6);

end
