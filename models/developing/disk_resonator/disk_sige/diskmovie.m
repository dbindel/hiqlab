% HiQLab
% Copyright (c): Regents of the University of California
% $Id: diskmovie.m,v 1.2 2006/05/04 05:37:44 tkoyama Exp $
clear;
qcloseall;

% -- Set parameters
param = [];
param.order = 3;
param.dense = 1.5e-6/4;
param.hdisk = 1.4e-6;

% -- Forcing frequency
w0 = 2*pi*47.5e6;

% -- Plotting parameters
plotopt = [];
plotopt.nframes = 60;
plotopt.fpcycle = 30;
plotopt.deform = 1e-6;
plotopt.cscale = 0;
plotopt.axis = [0, 50e-6, -5e-6, 5e-6];
plotopt.xlabel = {'', 'r (\mum)'};
plotopt.ylabel = 'z (\mum)';
plotopt.xscale = 1e6;
plotopt.titles = {'Radial displacement', 'Vertical displacement'};

% -- Load mesh from Lua input file
mesh  = Mesh_load('diskmesh.lua', param);

% -- Construct mass, stiffness matrix,
%    forcing vector and solve for
%    harmonic displacements.
drive_pat = Mesh_get_drive_f(mesh,'force_pattern');
harmonic_state(mesh,drive_pat,w0);

% -- Plot forced response
figure(1);
plotcycle2d(mesh,plotopt.deform,plotopt);
