% driver_force(wforce, params, plotopt)
%
% Driver for simulation of an excited Michigan disk resonator.
% Display time-harmonic behavior when forced at the specified freq.

function driver_force(wforce, params, plotopt)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_force.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $

if nargin < 2, params  = []; end
if nargin < 3, plotopt = []; end

mesh  = Mesh_load('diskmesh.lua', params);
[M,K] = Mesh_assemble_mk(mesh);
F     = Mesh_assemble_R(mesh);
Mesh_set_u(mesh, (K - wforce^2*M) \ F);

plotcycle2d(mesh,0,plotopt);
Mesh_delete(mesh);
