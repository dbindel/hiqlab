% plotfield1d(mesh, opt)
%
% Colors are chosen according to the magnitudes of the u components.
%
% The following options can be set via opt:
%   ufields  ([1  ]) - fields to use for displacement?
%   cfields  ([1  ]) - fields to use for color?
%   ncfields (1)     - number of color fields to obtain
%   axequal  (0)     - make axes equal?
%   subplot  ([])    - subplot size (default is [nfields 1])
%   deform   (0)     - Magnitude of deformation
%   clf      (1)     - Clear figure before plotting
%   xscale   (1)     - Amount to scale by (for unit change)
%   xlabel   ([])    - X label
%   ylabel   ([])    - Y label
%   titles   ([])    - Title (can be a cell array of titles)

function plotfield1d(mesh, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, opt = []; end

ufields  = qoptdefault(opt, 'ufields',  [1  ]);
cfields  = qoptdefault(opt, 'cfields',  [1  ]);
ncfields = qoptdefault(opt, 'ncfields', 1    );
axequal  = qoptdefault(opt, 'axequal',  0    );
subplm   = qoptdefault(opt, 'subplot',  [ncfields 1]);
deform   = qoptdefault(opt, 'deform',   0    );
clflag   = qoptdefault(opt, 'clf',      1    );
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
e_nen = zeros(1,numelt);
for i = 1:numelt
  e_nen(i) = Mesh_get_nen_elt(mesh,i);
end;
eplotf = e(:  ,find(e_nen==nen)); % Face  elements
eplotl = e(1:2,find(e_nen==2  )); % Line  elements
eplotp = e(1  ,find(e_nen==1  )); % Point elements
eplotl = eplotf;

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

% -- Evaluate function if character string
if ischar(cfields)
  uc = Mesh_get_lua_fields(mesh, cfields, ncfields)';
elseif ncfields > 0
  uc = real(u(cfields,:)).';
end

if clflag,                    clf;         end
if length(ufields) == 2,      pdisp = p + deform*real(u(ufields,:));
else,                         pdisp = p;
end

if ncfields > 0
  for j = 1:ncfields
    subplot(subplm(1),subplm(2),j);
    plot( pdisp,uc(:,j) , 'r*');
    if ischar(titles), title(titles);    end
    if iscell(titles), title(titles{j}); end
    if ischar(xlab),   xlabel(xlab);     end
    if iscell(xlab),   xlabel(xlab{j});  end
    if ischar(ylab),   ylabel(ylab);     end
    if iscell(ylab),   xlabel(ylab{j});  end
    if axequal,        axis equal;       end
  end
else
  plot( pdisp,uc(:,1) , 'r*');
  if ischar(titles), title(titles);    end
  if iscell(titles), title(titles{1}); end
  if ischar(xlab),   xlabel(xlab);     end
  if iscell(xlab),   xlabel(xlab{1});  end
  if ischar(ylab),   ylabel(ylab);     end
  if iscell(ylab),   xlabel(ylab{1});  end
  if axequal,        axis equal;       end
end

