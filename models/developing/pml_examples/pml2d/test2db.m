
% HiQLab
% Copyright (c): Regents of the University of California
% $Id: test2db.m,v 1.1 2006/05/01 08:20:43 tkoyama Exp $

try

  % -- Set up system parameters

  w = 1;
  xsym = 0;
  a = 1;
  b = 0;
  % -- Assemble and solve linear system

  mesh  = Mesh_load('test2db.lua');
  [M,K] = Mesh_assemble_mk(mesh);
  KK = K - w^2*M;

  p  = Mesh_get_x(mesh);
  e  = Mesh_get_e(mesh);
  id = Mesh_get_id(mesh);

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
  Mesh_set_u(mesh, u);

  % -- Plot forced response
  clear opt;
  opt.axis = [-12,12,-12,12];
  opt.nframes = 48;
  opt.npcycle = 24;
  plotcycle2d(mesh, 1, opt);
  Mesh_delete(mesh);

catch

  % -- Cleanup on error or break

  if (real(mesh) ~= 0)
    Mesh_delete(mesh);
  end
  error(lasterr);

end
