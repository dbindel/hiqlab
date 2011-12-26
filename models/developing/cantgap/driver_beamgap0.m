% HiQLab
% Copyright (c): Regents of the University of California
% $Id: driver_beamgap0.m,v 1.3 2005/10/31 12:12:57 tkoyama Exp $
 
param = [];
mesh = Mesh_load('beamgap.lua', param);
plotmesh(mesh);
axis equal;
pause;

u = zeros(Mesh_get_numid(mesh), 1);
Mesh_set_u(mesh,u);
for k = 1:5
  K  = Mesh_assemble_dR(mesh, 1, 0, 0);
  R  = Mesh_assemble_R(mesh);
  du = K\R;
  u  = u-du;
  Mesh_set_u(mesh, u);
  norm(R)
end

opt.cfields = [3];
opt.ufields = [1,2];
opt.deform  = 000;
opt.axequal = 1;
plotfield2d(mesh,opt);

l = 100;
g = 2;
w = 2;
x = Mesh_get_x(mesh);
nodeno = find(x(1,:) == l & x(2,:) == (g) );
disp = Mesh_get_disp(mesh);
nodeno
disp(:,nodeno)

Mesh_delete(mesh);
