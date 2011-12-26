% -- Parameters
param.f0    = 0;
param.order = 2;
param.dense = 5e-6;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('ring_only.lua',param);

% -- Display the mesh
figure(1);
opt.axequal = 1;
plotmesh(mesh,opt);

% -- Compute frequencies, Qs, and eigenvectors
%w0      = 618.50e6*2*pi;
%w0      = 705.50e6*2*pi;
w0      = 620.50e6*2*pi;
nev     = 25;
[V,w,Q] = mechmode(mesh, w0, nev);
freq = w/2/pi;

% -- Put obtained modes back into the mesh
for i = 1:nev
  show_mode = i;
  Mesh_set_u(mesh,V(:,show_mode));

  fprintf('Figure      :%d\n',i);
  fprintf('Freq [Hz]   :%d\n',real(freq(i)));
  fprintf('            :%d\n',imag(freq(i)));
  fprintf('Q           :%d\n',Q(i)         );

  % -- Plot mode shape
  figure(i);
  opt.deform  = 1e-5/Mesh_get_scale(mesh,'L');
  opt.axequal = 1;
  plotfield2d(mesh,opt);

  % -- Plot animation of modeshape
%  figure(3);
%  opt.cfields = [1,2];
%  opt.cbias   = 1;
%  plotcycle2d(mesh,opt.deform,opt);
end

