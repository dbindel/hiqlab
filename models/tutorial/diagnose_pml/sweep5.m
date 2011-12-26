% HiQLab
% Copyright (c): Regents of the University of California
% $Id: sweep5.m,v 1.1 2006/05/01 04:47:55 tkoyama Exp $

% Vary log(damping)  => dpw
% Vary number of elt
% Fix elt / wavelength

param.Ne1   = 10;

for i1 = 1:3
  
  figure(i1);
  clf;
  colormap(bone);
  hold on;
  
  for i2 = 1:4

    listN = 1:floor(48/i1);
    listD = 1:0.25:12;
    plotD = 1:0.5:8;
    
    param.order = i1;
    param.epw   = i2*(12/param.order);
    k           = 2*pi/param.epw;

    rl = [];
    for iN = 1:length(listN)
      for iD = 1:length(listD)
    
	Lpml = listN(iN);
	param.dpw = listD(iD) * log(10) * param.epw / k / Lpml^2;
	param.Ne2 = listN(iN);
	[rcoeff, kcoeff, resid] = diagnose_pml(param, k);
	rl(iN,iD) = rcoeff;
    
      end
    end
    
    subplot(2,2,i2);
    contour(listN*i1, listD, -log10(rl)', plotD);
    
  end
end

for i1 = 1:3
  
  figure(i1);
  subplot(2,2,1); ylabel('-Log10(Nominal reflection)');
  subplot(2,2,3); ylabel('-Log10(Nominal reflection)');
  subplot(2,2,3); xlabel('Nodes in PML');
  subplot(2,2,4); xlabel('Nodes in PML');
  
  for i2 = 1:4
    subplot(2,2,i2);
    param.epw = i2*(12/i1);
    title(sprintf('%d elts/wave (order %d)', param.epw, i1));
  end
  
  print('-deps', sprintf('sweep5_order%d.eps', i1));
  
end
