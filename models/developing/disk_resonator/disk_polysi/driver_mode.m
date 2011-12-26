% [w,Q,V] = driver_mode(w0, params)
%
% Compute the mode nearest shift w0.

function [w,Q,VV] = driver_mode(w0, params, nmode)
 
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_mode.m,v 1.1 2006/05/01 05:04:09 tkoyama Exp $

if nargin < 2, params = []; end
if nargin < 3, nmode = 1; end

mesh = Mesh_load('diskmesh.lua', params);
[M,K] = Mesh_assemble_mk(mesh);
Mesh_delete(mesh);

[VV,w,Q] = pml_mode(M,K,w0,nmode);

