% plotmesh(p, e, stretch)
%
% Plot the mesh described by the node position matrix p and element
% matrix e, with coloring according to the stretch function value.

function plotmesh_stretch(p, e, stretch)

% HiQLab
% Copyright (c): Regents of the University of California

if exist('OCTAVE_VERSION')
  disp('Octave does not have sufficiently MATLAB-compatible graphics');
  return;
end

figure(1);
eplot = plotelt2d(e);
col = sum(stretch,1);
caxis('auto');
patch('Vertices',   p', ...
      'Faces',      eplot', ...
      'FaceColor',  'interp', ...
      'FaceVertexCData', col');

