% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep1.m,v 1.1 2006/05/01 04:47:53 tkoyama Exp $

% Vary PML length and damping per wavelength
% Plot contours of damping

param.epw   = 10;
param.dpw   = 40;
param.Ne1   = param.epw;
param.Ne2   = 20;
param.order = 1;
k           = 2*pi/param.epw;

rl = [];
dpwl = [0.1:0.1:0.9, 1:20];
Ne2l = 1:60;

for idpw = 1:length(dpwl)
  for iNe2 = 1:length(Ne2l)
    param.dpw = dpwl(idpw);
    param.Ne2 = Ne2l(iNe2);
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(idpw,iNe2) = rcoeff;
  end
end

figure(1);
[cs,h] = contour(Ne2l, dpwl, log10(rl)); clabel(cs,h)
xlabel('Number of elements');
ylabel('Damping per wavelength');

figure(2);
semilogy(dpwl, rl(:,end));
xlabel('Damping per wavelength');
ylabel('Reflection coefficient');

