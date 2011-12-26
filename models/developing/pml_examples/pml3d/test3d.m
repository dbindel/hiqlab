% Driver for simple PML test: post sitting on a half space
% Axisymmetric version.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test3d.m,v 1.1 2006/05/01 08:20:56 tkoyama Exp $
 
rho = 1;
E   = 10;
nu  = 0.3;

k = 2;
a = 1;  % Vertical displacement of block
b = 0;  % Tilt of block
s = 1;  % Scale for plot displacement

f0 = 40;
m  = 1;

% --

try
  
  order = 2;
  nen   = (order+1)^3;

  mesh = Mesh_new(3);
  D = PMLElastic3d_new(E, nu, rho);
  Mesh_own_elt(mesh, D);

  Mesh_add_block3d(mesh, -10,-10,-10, 10,10,0, 11,11,5, D, order);
  p = Mesh_get_x(mesh);

  stretch = (abs(p)-6)*f0/k;
  stretch = stretch .* (stretch > 0);
  PMLElement_set_stretch(D, stretch);

  Mesh_initialize(mesh);

  [M,K] = Mesh_assemble_mk(mesh);
  KK    = K-k^2*M;
  id    = Mesh_get_id(mesh);

  Mesh_delete(mesh);
  
catch
  
  if real(mesh) ~= 0
    Mesh_delete(mesh);
  end
  error(lasterr);
  
end

% --

x = linspace(-10,10, 11);
y = linspace(-10,10, 11);
z = linspace(-10, 0,  5);

Jrigid  = find( sum(abs(p)) < 1e-8 );
Jrest   = find( sum(abs(p)) > 1e-8 );
Iactive = reshape(id(:,Jrest), 3*length(Jrest), 1);

xbc       = p(1,Jrigid);
ubc       = a + b*xbc;
Fbc       = -full(KK(:,id(3,Jrigid))*ubc);

u               = zeros(length(K),1);
u(id(3,Jrigid)) = ubc;
u(Iactive)      = KK(Iactive,Iactive) \ Fbc(Iactive);
u = reshape(u, 3, length(u)/3);

figure(1);
set(gcf, 'DoubleBuffer', 'on');

Iwant = find(p(2,:)==0);
for i = 1:80

  c = exp(complex(0, i*pi/12));
  i = i+1;

  surf( x, z, reshape( real(c*u(3,Iwant)),length(z),length(x) ) );
  axis([x(1), x(end), z(1), z(end), -1, 1]);
  caxis([-0.2, 0.2]);

  pause(0.25);
  
end

