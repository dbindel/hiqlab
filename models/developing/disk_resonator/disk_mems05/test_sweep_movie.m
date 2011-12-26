 
% Test case: Movie of frequency response vs. thickness

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_sweep_movie.m,v 1.1 2006/05/01 05:06:55 tkoyama Exp $

param = [];
param.order = 3;
param.dense = 1e-6/4;
w0 = 2*pi*47.5e6;
t = linspace(1.4e-6, 1.6e-6, 21);

plotopt.framepng = '/tmp/image%4.4d';
plotopt.nframes = 50;
plotopt.fpcycle = 50;
plotopt.cscale = 1;
[w,Q] = driver_sweep(w0, 2, 'hdisk', t, param, plotopt);
