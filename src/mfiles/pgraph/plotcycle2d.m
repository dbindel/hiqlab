% plotcycle2d(mesh, s, opt)
%
% Plot x and y components for a time-harmonic response.  Colors are
% chosen according to the magnitudes of the u components, and the actual
% displayed wave amplitude is scaled by s.  If it is not provided,
% s defaults to 1.
%
% The following options can be set via opt:
%   framepng ([])    - if not empty, a frame file format string
%   nframes  (32)    - number of frames to be plotted
%   startf   (1)     - start frame number
%   fpcycle  (16)    - frames per cycle
%   fpause   (0.25)  - pause time betseen frames
%   cscale   (0)     - scale x and y colors together?
%   cbias    (3)     - bias of the color scale (cmax/cbias->red)
%   ufields  ([1 2]) - fields to use for displacement?
%   cfields  ([1 2]) - fields to use for color?
%   axequal  (0)     - make axes equal?
%   axis     ([])    - set axes
%   subplot  ([])    - subplot size (default is [nfields 1])
%   xscale   (1)     - Amount to scale by (for unit change)
%   xlabel   ([])    - X label
%   ylabel   ([])    - Y label
%   titles   ([])    - Title (can be a cell array of titles)
%   avi_file ([])    - Title of avi file to generate
%                      (Default: NO AVI FILE IS PRODUCED)
%   avi_w    (300)   - Width  of the window size
%   avi_h    (300)   - Height of the window size
%   avi_left ( 10)   - Window distance from left side of monitor
%   avi_bottom(def)  - Window distance from bottom    of monitor
%                   (Default is set so window touches top of monitor)
function plotcycle2d(mesh, s, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, s   = 1;  end
if nargin < 3, opt = []; end

framepng = qoptdefault(opt, 'framepng', []   );
nframes  = qoptdefault(opt, 'nframes',  32   );
sframe   = qoptdefault(opt, 'sframe',   1    );
fpcycle  = qoptdefault(opt, 'fpcycle',  16   );
fpause   = qoptdefault(opt, 'fpause',   0.25 );
cscale   = qoptdefault(opt, 'cscale',   0    );
cbias    = qoptdefault(opt, 'cbias',    3    );
ufields  = qoptdefault(opt, 'ufields',  [1 2]);
cfields  = qoptdefault(opt, 'cfields',  [1 2]);
axequal  = qoptdefault(opt, 'axequal',  0    );
ax       = qoptdefault(opt, 'axis',     []   );
subplm   = qoptdefault(opt, 'subplot',  [length(cfields) 1]);
xscale   = qoptdefault(opt, 'xscale',   1    );
xlab     = qoptdefault(opt, 'xlabel',   []   );
ylab     = qoptdefault(opt, 'ylabel',   []   );
titles   = qoptdefault(opt, 'titles',   []   );

ax = ax*abs(xscale);
p = Mesh_get_x(mesh) * abs(xscale);
e = Mesh_get_e(mesh);
u = Mesh_get_disp(mesh) * xscale;

% -- Pick only those elements which have nen nodes per element
numelt= Mesh_numelt(mesh);
nen   = Mesh_get_nen(mesh);
e_nen = Mesh_get_nen_elt(mesh);
eplotf = e(:  ,find(e_nen==nen)); % Face  elements
eplotl = e(1:2,find(e_nen==2  )); % Line  elements
eplotp = e(1  ,find(e_nen==1  )); % Point elements
eplotf = plotelt2d(eplotf);

if length(cfields) > 0
  cmax = max(abs(u(cfields,:)), [], 2);
  if cscale
    cmax = ones(size(cmax)) * max(cmax);
  end
end

boxx = [];
pu = [];

if length(ufields) == 2
  boxx = [ min(p(1,:)-s*abs(u(ufields(1),:))), ...
           min(p(2,:)-s*abs(u(ufields(2),:))), ...
           max(p(1,:)+s*abs(u(ufields(1),:))), ...
           max(p(2,:)+s*abs(u(ufields(2),:))) ];
  pu = s*u(ufields,:);
end
 
if length(cfields) > 0
  pgraph_open(length(cfields));
  for j = 1:length(cfields)
    pgraph_current(j);
    if ~isempty(boxx)
      pgraph_bound(boxx(1), boxx(2), boxx(3), boxx(4));
    end
    pgraph_write(p, eplotf);
    pgraph_write_mode(p, pu, u(cfields(j),:));
    pgraph_cbound(-cmax(j), cmax(j));
    pgraph_animate(1);
  end
else
  pgraph_open();
  pgraph_current(j);
  if ~isempty(boxx)
    pgraph_bound(boxx(1), boxx(2), boxx(3), boxx(4));
  end
  pgraph_write(p, eplotf);
  pgraph_write_mode(p, pu);
  pgraph_animate(1);
end

