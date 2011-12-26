% plotfield2d(mesh, opt)
%
% Colors are chosen according to the magnitudes of the u components.
%
% The following options can be set via opt:
%   cscale   (0)     - scale field colors together?
%   cbias    (3)     - bias of the color scale (cmax/cbias->red)
%   ufields  ([1 2]) - fields to use for displacement?
%   cfields  ([1 2]) - fields to use for color?
%   ncfields (1)     - number of color fields to obtain
%   axequal  (0)     - make axes equal?
%   subplot  ([])    - subplot size (default is [nfields 1])
%   deform   (0)     - Magnitude of deformation
%   clf      (1)     - Clear figure before plotting
%   cbar     (0)     - Add colorbar
%   xscale   (1)     - Amount to scale by (for unit change)
%   xlabel   ([])    - X label
%   ylabel   ([])    - Y label
%   titles   ([])    - Title (can be a cell array of titles)

function plotfield2d(mesh, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if exist('OCTAVE_VERSION')
  disp('Octave does not have sufficiently MATLAB-compatible graphics');
  return;
end

if nargin < 2, opt = []; end

cscale   = qoptdefault(opt, 'cscale',   0    );
cbias    = qoptdefault(opt, 'cbias',    3    );
ufields  = qoptdefault(opt, 'ufields',  [1 2]);
cfields  = qoptdefault(opt, 'cfields',  [1 2]);
if ischar(cfields)
  ncfields = qoptdefault(opt, 'ncfields', 1);
else
  ncfields = qoptdefault(opt, 'ncfields', length(cfields));
end
axequal  = qoptdefault(opt, 'axequal',  0    );
subplm   = qoptdefault(opt, 'subplot',  [ncfields 1]);
deform   = qoptdefault(opt, 'deform',   0    );
clflag   = qoptdefault(opt, 'clf',      1    );
cbar     = qoptdefault(opt, 'cbar',     0    );
xscale   = qoptdefault(opt, 'xscale',   1    );
xlab     = qoptdefault(opt, 'xlabel',   []   );
ylab     = qoptdefault(opt, 'ylabel',   []   );
titles   = qoptdefault(opt, 'titles',   []   );

p = Mesh_get_x(mesh) * abs(xscale);
u = Mesh_get_disp(mesh) * xscale;
e = Mesh_get_e(mesh);

% -- Pick only those elements which have nen nodes per element
numelt= Mesh_numelt(mesh);
nen   = Mesh_get_nen(mesh);
e_nen = Mesh_get_nen_elt(mesh);
eplotf = e(:  ,find(e_nen==nen)); % Face  elements
eplotl = e(1:2,find(e_nen==2  )); % Line  elements
eplotp = e(1  ,find(e_nen==1  )); % Point elements
eplotf = plotelt2d(eplotf);

if ischar(cfields)
  uc = Mesh_get_lua_fields(mesh, cfields, ncfields)';
  cmax = max(abs(uc));
  if cscale
    cmax = ones(size(cmax)) * max(cmax);
  end
elseif ncfields > 0
  uc = real(u(cfields,:)).';
  cmax = max(abs(u(cfields,:)), [], 2);
  if cscale
    cmax = ones(size(cmax)) * max(cmax);
  end
end

if clflag,                    clf;         end
if length(ufields) == 2,      pdisp = p + deform*real(u(ufields,:));
else,                         pdisp = p;
end

if ncfields > 0
  for j = 1:ncfields
    subplot(subplm(1),subplm(2),j);
    patch('Vertices', pdisp.', ...
          'Faces',    eplotf', ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'interp', ...
          'FaceVertexCData', uc(:,j));
    caxis([-cmax(j), cmax(j)]/cbias);
    if ischar(titles), title(titles);    end
    if iscell(titles), title(titles{j}); end
    if ischar(xlab),   xlabel(xlab);     end
    if iscell(xlab),   xlabel(xlab{j});  end
    if ischar(ylab),   ylabel(ylab);     end
    if iscell(ylab),   xlabel(ylab{j});  end
    if axequal,        axis equal;       end
    if cbar,           colorbar;         end
  end
else
  patch('Vertices', pdisp.', ...
        'Faces',    eplotf', ...
        'FaceColor', 'b');
  if ischar(titles), title(titles);    end
  if iscell(titles), title(titles{1}); end
  if ischar(xlab),   xlabel(xlab);     end
  if iscell(xlab),   xlabel(xlab{1});  end
  if ischar(ylab),   ylabel(ylab);     end
  if iscell(ylab),   xlabel(ylab{1});  end
  if axequal,        axis equal;       end
end
