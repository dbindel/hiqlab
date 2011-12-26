% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep6.m,v 1.1 2006/05/01 04:47:49 tkoyama Exp $

% Vary damping per element and number of elements;
% compare to max of nominal and long-range damping

param.epw   = 30;
param.dpw   = 40;
param.Ne1   = 10;
param.Ne2   = 100;
param.order = 3;
k           = 2*pi/param.epw;

rl   = [];
rnl  = [];
dpel = logspace(-2,2,100);
nel  = 1:30;

for idpe = 1:length(dpel)
  for ine = 1:length(nel)
    param.dpw = dpel(idpe) * param.epw;
    param.Ne2 = nel(ine);
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpe,ine)  = rcoeff;
    rnl(idpe,ine) = exp(-dpel(idpe) * k * nel(ine)^2);
  end
end

figure(1);
[cs,h] = contour(nel, log10(dpel), log10(rl)); 
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
title('Computed damping profile');
clabel(cs,h);

figure(2);
rel = rl(:,end)*ones(1,length(nel));
[cs,h] = contour(nel, log10(dpel), max(log10(rnl),log10(rel)) ); 
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
title('Predicted damping profile');
clabel(cs,h);

figure(3);
dd = log10(rl) - log10(rnl+rel)/;
[cs,h] = contour(nel, log10(dpel), dd);
hold on
[cs,h] = contour(nel, log10(dpel), log10(rnl)-log10(rel), [0 0]);
set(h, 'edgecolor', 'b');
set(h, 'linestyle', ':');
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
title(sprintf('Difference from predicted damping: %f to %f', ...
	      min(min(dd)), max(max(dd)) ));

figure(4);
surf(nel, log10(dpel), log10(rl));
xlabel('Number of elements');
ylabel('Log10(Damping per element)');
zlabel('Log10(Reflection coeff)');
