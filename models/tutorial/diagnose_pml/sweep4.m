% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep4.m,v 1.1 2006/05/01 04:47:48 tkoyama Exp $

% Vary log(damping)  => dpw
% Vary number of elt
% Fix elt / wavelength

param.epw   = 8;
param.order = 3;
param.Ne1   = 10;
k           = 2*pi/param.epw;

rl = [];
listN = 1:50;
listD = 1:0.5:12;


for iN = 1:length(listN)
  for iD = 1:length(listD)
    
    Lpml = listN(iN);
    param.dpw = listD(iD) * log(10) * param.epw / k / Lpml^2;
    param.Ne2 = listN(iN);
    [rcoeff, kcoeff, resid] = diagnose_pml(param, k);
    rl(iN,iD) = rcoeff;
    
  end
end

figure(1);
[cs,h] = contour(listN, listD, log10(rl)'); 
xlabel('Number of elements');
ylabel('Log10(Nominal damping)');
%clabel(cs,h);

figure(2);
surf(listN, listD, log10(rl)');
