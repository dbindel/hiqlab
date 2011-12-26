% plotmesh1d(mesh, opt)
%
% Plot the 1d mesh.
%
% The following options can be set via opt:
%   anchors ('g+') - Mark nodes with displacement BC (empty for no marker)
%   forces  ('r+') - Mark nodes with force BC (empty for no marker)
%   deform  ( 0  ) - Deform mesh?
%   clf     ( 1  ) - Clear the figure before display?
%   axequal ( 0  ) - make axes equal?

function plotmesh1d(mesh, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, opt = []; end

anchors = qoptdefault(opt, 'anchors', 'g+');
forces  = qoptdefault(opt, 'forces',  'r+');
deform  = qoptdefault(opt, 'deform',   0  );
do_clf  = qoptdefault(opt, 'clf',      1  );
axequal = qoptdefault(opt, 'axequal',  0  );

p = Mesh_get_x(mesh);
e = Mesh_get_e(mesh);

% -- Pick only those elements which have nen nodes per element
numelt= Mesh_numelt(mesh);
nen   = Mesh_get_nen(mesh);
e_nen = zeros(1,numelt);
for i = 1:numelt
  e_nen(i) = Mesh_get_nen_elt(mesh,i);
end;
eplotf = e(:  ,find(e_nen==nen)); % Face  elements
eplotl = e(1:2,find(e_nen==2  )); % Line  elements
eplotp = e(1  ,find(e_nen==1  )); % Point elements
eplotl = eplotf;

if deform,
  u = Mesh_get_disp(mesh);
  u = u(1,:);
  p = p + u;
end
if do_clf,  clf;  end

% -- Plot nodes
plot( p, zeros(length(p)), 'ks')

% -- Plot Line elements
if ~isempty(eplotl)
  hold_save = ishold;
  hold on;
  for i = 1:size(eplotl,2)
    xy = p(:, eplotl(:,i)');
    plot(xy, zeros(length(xy)), 'k-')
  end;
  if ~hold_save, hold off; end;
end

if ~isempty(anchors) | ~isempty(forces)

  hold_save = ishold;
  hold on

  bc = Mesh_get_bc(mesh);
  [ndf,numnp] = size(bc);
  I = [];
  J = [];
  for i = 1:ndf
    I = union(I,find(bc(i,:)==1));
    J = union(J,find(bc(i,:)==2));
  end;
  if ~isempty(anchors), plot(p(1,I), zeros(length(I)), anchors); end
  if ~isempty(forces),  plot(p(1,J), zeros(length(J)), forces);  end

  if ~hold_save, hold off; end

end

if axequal,        axis equal;       end
