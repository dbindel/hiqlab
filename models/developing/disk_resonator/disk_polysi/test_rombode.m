 
% Test case: Run full model and reduced-order model for an example
% problem (diamond disk resonator).  Then visually compare.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_rombode.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $
 
% -- Set parameters

param.order = 3;
param.dense = 1e-6/8;
param.disk_material = 'diamond';
param.rpost =  1.6e-6 / 2;        % Post radius
param.rdisk = 24.0e-6 / 2;        % Disk radius
param.hdisk =  2.5e-6;            % Thickness of disk layer
param.hpost =  8e-7;              % Thickness of post layer

kmax = 3;         % Order of reduced model
npts = 1001;      % Number of points on Bode plot (reduced model)
Npts = 51;        % Number of points on Bode plot (full model)
freqlo = 453e6;   % Min frequency (Hz)
freqhi = 454e6;   % Max frequency (Hz)
w0     = 455e6;   % Shift frequency (Hz)

% -- Compute Bode plots

freq  = 2*pi*linspace(freqlo, freqhi, npts);
freqN = 2*pi*linspace(freqlo, freqhi, Npts);
w0    = 2*pi*w0;

tic; Hn = driver_bode(freq,  param, kmax, w0); toc
tic; HN = driver_bode(freqN, param);           toc

% -- Draw Bode plots

figure(1); clf; hold on
opt.lstyle = 'b';   plot_bode(freqN, HN, opt);
opt.lstyle = 'r--'; plot_bode(freq,  Hn, opt);
hold off
