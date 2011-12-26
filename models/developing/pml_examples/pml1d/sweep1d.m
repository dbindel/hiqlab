k = 2*pi/10;
dpw = linspace(10,500,50);

for j = 1:length(dpw)
  param.dpw = dpw(j);
  mesh  = Mesh_load('pml1d.lua', param);
  Mesh_make_harmonic(mesh, k);
  [M,K] = Mesh_assemble_mk(mesh);
  F     = Mesh_assemble_R(mesh);
  u     = -(K-k^2*M)\F;
  Mesh_set_u(mesh, u);

  x = Mesh_get_x(mesh);
  u = Mesh_get_disp(mesh);
  plot(x, real(u), x, imag(u)); axis tight
  title(sprintf('Damping param = %f', param.dpw));
  pause(0.4);

  Mesh_delete(mesh);
end
