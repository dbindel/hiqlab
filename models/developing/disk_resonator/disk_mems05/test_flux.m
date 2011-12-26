% driver_force(wforce, params, plotopt)
%
% Driver for simulation of an excited Michigan disk resonator.
% Display time-harmonic behavior when forced at the specified freq.

%function test_flux(wforce, params, plotopt)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test_flux.m,v 1.1 2006/05/01 05:06:55 tkoyama Exp $

%if nargin < 2, params  = []; end
%if nargin < 3, plotopt = []; end

params = [];
plotopt = [];

params.hdisk = 1.6e-6;
params.rdisk = 40e-6;
params.dense = 1e-6/2;
params.order = 3;
wforce = 45.04 * 2e6*pi;
plotopt.xscale = 1e-6;
plotopt.yscale = 1e-6;

mesh  = Mesh_load('diskmesh.lua', params);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
Mesh_set_u(mesh, -(K-wforce^2*M) \ F);
Mesh_make_harmonic(mesh, wforce);

p = Mesh_get_x(mesh);
E = Mesh_mean_power(mesh);

%I = 1:8:size(p,2);
%E = (2*pi*[1;1]*p(1,:)) .* E;
I = 1:size(p,2);
h = quiver(p(1,I), p(2,I), E(1,I), E(2,I));
axis equal;

% plotcycle2d(mesh,0,plotopt);
% Mesh_delete(mesh);
