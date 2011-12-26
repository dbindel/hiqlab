% plotmesh(mesh, opt)
%
% Plot the mesh.
%
% The following options can be set via opt:
%   anchors ('g+') - Mark nodes with displacement BC (empty for no marker)
%   forces  ('r+') - Mark nodes with force BC (empty for no marker)
%   deform  ( 0  ) - Deform mesh?
%   clf     ( 1  ) - Clear the figure before display?
%   axequal ( 0  ) - make axes equal?

function plotmesh(mesh, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, opt = []; end

ndm = Mesh_get_ndm(mesh);

if ndm == 1
   plotmesh1d(mesh,opt);
elseif ndm == 2
   plotmesh2d(mesh,opt);
else
   error('Plotmesh not supported for this mesh dimension');
end
