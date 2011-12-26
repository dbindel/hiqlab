clear;
qcloseall;
% -- Parameters
param.use_case  = 'using_electrodes_dyn.lua';
param.order     = 1;
param.dense     = 1.0e-6;

param.fill_gap  = 1;
param.reps      = 160;
param.Vf_dc     = 100;

param.f0        =   0;
f0              =  20;

param.mech      = 0;
nev             = 2;
w0              = 140e6*2*pi;

% -- Load mesh from Lua input file
[mesh, L] = Mesh_load('dielectric_drive.lua',param);

% -- Compute for static state (f0 = 0)
sopt.nonlinear = 'NR';
static_state(mesh,sopt);

% -- Check for pull-in
U              = Mesh_get_u(mesh);
sense_gap_dist = Mesh_get_sense_u(mesh,'sense_gap_dist');
Dt             = Lua_get_double(L,'Dt');
fprintf('Gap distance / 3     [um]:%d\n',               Dt/3/1e-6);
fprintf('Gap distance decrease[um]:%d\n',(sense_gap_dist'*U)/1e-6);

% -- Reload mesh with PML
Mesh_delete(mesh);
param.f0  = f0;
[mesh, L] = Mesh_load('dielectric_drive.lua',param);
Mesh_set_u(mesh,U);

% -- Extract element & global numbers
eno_d     = Lua_get_double(L,'eno_d')+1;
eno_s     = Lua_get_double(L,'eno_s')+1;
idg_d     = Lua_get_double(L,'idg_d')+1;
idg_s     = Lua_get_double(L,'idg_s')+1;
param.eno = [eno_d,eno_s];
param.idg = [idg_d,idg_s];

% -- Compute eigenvalues
[V,w,Q] = emcmode(mesh, w0, nev, param);
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
