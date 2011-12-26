% plot_mode(mesh,V,w,opt)
%
% Show mode shapes given in V in order. Optionally
% the animation can also be viewed.
%
% The options that can be set in opt follow those of
%
% - plotfield
% - plotcycle
%
% See these function for options.
%
% The following options can be set via opt:
%   animate(0) - Show animated plot ?
%
function plot_mode(mesh,V,w,opt)

if nargin < 4, opt = []; end

animate = qoptdefault(opt, 'animate',   0    );
Q       = abs(w)./(2*imag(w));
freq    = w/2/pi;
ndm     = Mesh_get_ndm(mesh);
nev     = length(freq);

% -- Put obtained modes back into the mesh
for i = 1:nev
  show_mode = i;
  Mesh_set_u(mesh,V(:,show_mode));

  fprintf('Figure      :%d\n',i);
  fprintf('Freq [Hz]   :%d\n',real(freq(i)));
  fprintf('            :%d\n',imag(freq(i)));
  fprintf('Q           :%d\n',Q(i)         );

  if ~animate
    % -- Plot mode shape
    if ndm==2,  plotfield2d(mesh,opt); end;
  else
    % -- Plot animation of modeshape
    if ndm==2,  plotcycle2d(mesh,opt.deform,opt); end;
  end
  if i~=nev, pause; end;
end

