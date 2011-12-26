% Driver for simple PML test: block sitting on a half space.
% Version 1: Don't use Lua for anything.

% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test2d.m,v 1.1 2006/05/01 08:20:48 tkoyama Exp $

try
 
  % -- Set up properties

  rho = 1;
  E   = 10;
  nu  = 0.3;

  xrad = 12;  % Radius of domain (x dimension)
  yrad = 12;  % Radius of domain (y dimension)
  xpml =  9;  % Start of PML (x dimension)
  ypml =  9;  % Start of PML (y dimension)
  xsym =  0;  % Symmetry boundary at x = 0?
  mx   = 36;  % Element spaces in x
  my   = 24;  % Element spaces in y
  f0   = 40;  % Maximum stretch parameter

  w = 1;      % Drive frequency
  a = 1;      % Magnitude of vertical pull in forcing
  b = 0;      % Magnitude of tilt in forcing
  s = 1;      % Stretch magnitude for displacement plot

  order = 3;           % Element order
  nen   = (order+1)^2; % Number of element nodes 

  % -- Assemble mesh quantities

  mesh = Mesh_new(2,nen,2);
  D    = PMLElastic2d_new(E, nu, rho, 0, nen);
  Mesh_own_elt(mesh, D);

  if xsym
    Mesh_add_block2d(mesh,     0,-yrad, xrad,0, mx+1,my+1, D, order);
  else
    Mesh_add_block2d(mesh, -xrad,-yrad, xrad,0, mx+1,my+1, D, order);
  end
  numnp  = Mesh_numnp(mesh);
  numelt = Mesh_numelt(mesh);

  p = Mesh_get_x(mesh);
  e = Mesh_get_e(mesh)+1;

  stretch = [(abs(p(1,:))-xpml)/(xrad-xpml); ...
	     (abs(p(2,:))-ypml)/(yrad-ypml)] * f0 / w;
  stretch = stretch .* (stretch > 0);
  PMLElement_set_stretch(D, stretch);

  Mesh_initialize(mesh);

  [M,K] = Mesh_assemble_mk(mesh);
  id    = Mesh_get_id(mesh)+1;

  KK    = K - w^2*M;

  % -- Find where to apply BC

  Jrigid  = find( p(2,:) == 0 & abs(p(1,:)) <= 1 );
  Jzerox  = find( p(1,:) == 0 );

  idred = id;
  idred(:,Jrigid) = 0;
  if xsym, idred(1,Jzerox) = 0; end
  Iactive = idred(find(idred));

  % -- Compute boundary forcing and solve

  xbc = p(1,Jrigid)';
  ubc = a + b*xbc;
  Fbc = -KK(:,id(2,Jrigid))*ubc;

  u               = zeros(length(K),1);
  u(id(2,Jrigid)) = ubc;
  u(Iactive)      = KK(Iactive,Iactive) \ Fbc(Iactive);

  % --- Plot forced response

  Mesh_set_u(mesh, u);
  plotcycle2d(mesh, s);
  Mesh_delete(mesh);

catch

  if real(mesh) ~= 0
    Mesh_delete(mesh);
  end
  error(lasterr);

end

