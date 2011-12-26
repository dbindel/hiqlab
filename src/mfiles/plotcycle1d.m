% plotcycle1d(mesh, s, opt)
%
% Plot a time-harmonic response.  The actual displayed wave
% amplitude is scaled by s.  If it is not provided, s defaults to 1.
%
% The following options can be set via opt:
%   framepng ([])    - if not empty, a frame file format string
%   nframes  (32)    - number of frames to be plotted
%   startf   (1)     - start frame number
%   fpcycle  (16)    - frames per cycle
%   fpause   (0.25)  - pause time betseen frames
%   axis     ([])    - set axes
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

function plotcycle1d(mesh, s, opt)

% HiQLab
% Copyright (c): Regents of the University of California

if nargin < 2, s   = 1;  end
if nargin < 3, opt = []; end

framepng = qoptdefault(opt, 'framepng', []   );
nframes  = qoptdefault(opt, 'nframes',  32   );
sframe   = qoptdefault(opt, 'sframe',   1    );
fpcycle  = qoptdefault(opt, 'fpcycle',  16   );
fpause   = qoptdefault(opt, 'fpause',   0.25 );
ax       = qoptdefault(opt, 'axis',     []   );
xscale   = qoptdefault(opt, 'xscale',   1    );
xlab     = qoptdefault(opt, 'xlabel',   []   );
ylab     = qoptdefault(opt, 'ylabel',   []   );
titles   = qoptdefault(opt, 'titles',   []   );
avi_file = qoptdefault(opt, 'avi_file', []   );
avi_w    = qoptdefault(opt, 'avi_w'   , 300  );
avi_h    = qoptdefault(opt, 'avi_h'   , 300  );
monitor_size = get( 0, 'ScreenSize');
avi_left  = qoptdefault(opt, 'avi_left'  ,  10  );
avi_bottom= qoptdefault(opt, 'avi_bottom', monitor_size(4)-avi_h );

if ~isempty(avi_file)
  fig = gcf;
  set(fig,'DoubleBuffer','on');
  set(gca,'Visible','off');
  set(gcf,'Position',[avi_left, avi_bottom, avi_w, avi_h]);
  mov = avifile(avi_file);
end

ax = ax*abs(xscale);
p = Mesh_get_x(mesh) * abs(xscale);
e = Mesh_get_e(mesh);
u = Mesh_get_disp(mesh) * xscale;
u = u(1,:);

for i = 1:nframes

  clf;
  c = exp(complex(0, i*2*pi/fpcycle));

  line(p(e), s*real(c*u(e)), 'color', 'b');

  if ischar(titles), title(titles);    end
  if iscell(titles), title(titles{j}); end
  if ischar(xlab),   xlabel(xlab);     end
  if iscell(xlab),   xlabel(xlab{j});  end
  if ischar(ylab),   ylabel(ylab);     end
  if iscell(ylab),   xlabel(ylab{j});  end
  if ~isempty(ax), axis(ax);   end

  if ~isempty(avi_file)
    F   = getframe(gcf);
    mov = addframe(mov,F);
    pause(fpause*2);
  end

  if ~isempty(framepng)
    print('-dpng', '-r0', sprintf(framepng, sframe + i-1));
  else
    pause(fpause);
  end

end

if ~isempty(avi_file)
  mov = close(mov);
end
