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
%   xscale   (1)   - Amount to scale by (for unit change)

function plotmesh2d(mesh, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, opt = []; end

anchors = qoptdefault(opt, 'anchors', 'g+');
forces  = qoptdefault(opt, 'forces',  'r+');
deform  = qoptdefault(opt, 'deform',   0  );
do_clf  = qoptdefault(opt, 'clf',      1  );
axequal = qoptdefault(opt, 'axequal',  0  );
xscale  = qoptdefault(opt, 'xscale',   1  );

p = Mesh_get_x(mesh) * abs(xscale);
e = Mesh_get_e(mesh);

% -- Pick only those elements which have nen nodes per element
numelt= Mesh_numelt(mesh);
nen   = Mesh_get_nen(mesh);
e_nen = Mesh_get_nen_elt(mesh);
eplotf = e(:  ,find(e_nen==nen)); % Face  elements
eplotl = e(1:2,find(e_nen==2  )); % Line  elements
eplotp = e(1  ,find(e_nen==1  )); % Point elements

% -- Plot Face elements
if deform,
  u = Mesh_get_disp(mesh) * xscale;
  u = u(1:2,:);
  p = p + u;
  p = real(p);
end
% if do_clf,  clf;  end
if ~isempty(eplotf) & (nen>=4)
  eplotf = plotelt2d(eplotf);
  pgraph_setskeleton(1);
  pgraph_write(p, eplotf);
end;

% --------------------------------------------
% FIXME: These things aren't implemented yet

if 0

% -- Plot Line elements
if ~isempty(eplotl)
  hold_save = ishold;
  hold on;
  for i = 1:size(eplotl,2)
    xy = p(:, eplotl(:,i)');
    plot(xy(1,:), xy(2,:), 'rs-',...
                           'LineWidth',4,...
                           'MarkerEdgeColor','r',...
                           'MarkerFaceColor','r',...
                           'MarkerSize',12);
  end;
  if ~hold_save, hold off; end;
end

% -- Plot Point elements
if ~isempty(eplotp)
  hold_save = ishold;
  hold on;
  for i = 1:size(eplotp,2)
    xy = p(:, eplotp(:,i)');
    plot(xy(1,:), xy(2,:), 'ms',...
                           'LineWidth',4,...
                           'MarkerEdgeColor','m',...
                           'MarkerFaceColor','m',...
                           'MarkerSize',24);
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
  if ~isempty(anchors), plot(p(1,I), p(2,I), anchors); end
  if ~isempty(forces),  plot(p(1,J), p(2,J), forces);  end

  if ~hold_save, hold off; end

end

if axequal,        axis equal;       end

end
