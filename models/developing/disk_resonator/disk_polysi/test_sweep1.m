 
% Compare sweeps over a variety of frequency ranges for different
% values of the PML parameters.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_sweep1.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
param.order = 3;
param.dense = 1e-6/6;
param.disk_material = 'diamond';
param.rpost =  1.6e-6 / 2;
param.rdisk = 24.0e-6 / 2;
param.hdisk =  3e-6;
param.hpost =  8e-7;
w0 = 2*pi*455e6;
t  = linspace(2e-6, 3e-6, 51);
  
param.rbd  = 2e-6;
param.rpml = 4e-6;
param.f0   = 20;
[w1, Q1] = driver_sweep(w0, 2, 'hdisk', t, param);

param.rbd  = 3e-6;
param.rpml = 5e-6;
param.f0   = 20;
[w2, Q2] = driver_sweep(w0, 2, 'hdisk', t, param);

param.rbd  = 3e-6;
param.rpml = 6e-6;
param.f0   = 20;
[w3, Q3] = driver_sweep(w0, 2, 'hdisk', t, param);

param.rbd  = 2e-6;
param.rpml = 4e-6;
param.f0   = 30;
[w4, Q4] = driver_sweep(w0, 2, 'hdisk', t, param);

param.rbd  = 2e-6;
param.rpml = 4e-6;
param.f0   = 40;
[w5, Q5] = driver_sweep(w0, 2, 'hdisk', t, param);

semilogy(t/1e-6, [max(Q1); max(Q2); max(Q3); max(Q4); max(Q5)]);
xlabel('Film thickness (um)');
ylabel('Q factor');
title('Frequency sweep for different PML settings');
