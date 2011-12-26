 
% Test case: Movie of frequency response vs. thickness

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_sweep2.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
param.order = 3;
param.dense = 1e-6/6;
param.disk_material = 'diamond';
param.rpost =  1.6e-6 / 2;
param.rdisk = 24.0e-6 / 2;
param.hdisk =  3e-6;
param.hpost =  8e-7;
w0 = 2*pi*455e6;
t = linspace(2e-6, 3e-6, 41);

plotopt.framepng = '/tmp/image%4.4d';
plotopt.nframes = 50;
plotopt.fpcycle = 50;
plotopt.cscale = 1;
[w,Q] = driver_sweep(w0, 2, 'hdisk', t, param, plotopt);
