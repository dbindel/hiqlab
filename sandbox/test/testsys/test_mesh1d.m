% Test mesh construction

% -- 1D linear field test
n = 10;
for order = 1:3

  % Construct a mesh
  mesh = Mesh_new(1);
  e = PMLScalar1d_new(123,45);
  Mesh_own_elt(mesh,e);
  Mesh_add_block1d(mesh,0,1,order*n+1, e,order);
  Mesh_add_block1d(mesh,1,3,order*n+1, e,order);
  Mesh_tie(mesh,1e-3);
  Mesh_initialize(mesh);

  % Construct a second mesh using add_nodes and add_elt; does it work?
  p1 = Mesh_get_x(mesh);
  e1 = Mesh_get_e(mesh);
  id1 = Mesh_get_id(mesh);
  bc1 = Mesh_get_bc(mesh);
  numnp1 = Mesh_get_numid(mesh);

  mesh2 = Mesh_new(1);
  Mesh_add_node(mesh2, p1);
  Mesh_add_element(mesh2, e1, e);
  Mesh_initialize(mesh2);
  p2 = Mesh_get_x(mesh);
  e2 = Mesh_get_e(mesh);
  id2 = Mesh_get_id(mesh);
  bc2 = Mesh_get_bc(mesh);
  numnp2 = Mesh_get_numid(mesh);
  Mesh_delete(mesh2);

  qassert(norm(p1-p2,1) == 0, 'Position vectors unequal!');
  qassert(norm(e1-e2,1) == 0, 'Element structures unequal!');
  qassert(norm(id1-id2,1) == 0, 'ID structures unequal!');
  qassert(numnp1 == numnp2, 'Numnp inconsistent!');

  % Check the scales (should all be one!)
  [cu,cf] = Mesh_get_scales(mesh);
  cL = Mesh_get_scale(mesh, 'L');
  qassert(cu == 1, 'Incorrect displacement scale');
  qassert(cf == 1, 'Incorrect force scale');
  qassert(cL == 1, 'Incorrect length scale');

  % Check the connectivity
  numelt = Mesh_numelt(mesh);
  nen    = Mesh_get_nen(mesh);
  qassert(numelt == 2*n, 'Wrong number of elements!');
  qassert(nen == order+1, 'Wrong number of element nodes!');
  lastx = -100;
  for j = 1:numelt
    nen = Mesh_get_nen_elt(mesh,j);
    qassert(nen == order+1, 'Wrong number of element nodes!');
    for i = 1:nen
      x = Mesh_x(mesh, 1, Mesh_ix(mesh,i,j));
      qassert(lastx <= x, 'Out-of-order nodes?');
      lastx = x;
    end     
  end

  % Get the mesh id assignment and test consistency
  id = Mesh_get_id(mesh);
  qassert(length(id) == Mesh_numnp(mesh), 'ID is the wrong length');

  % Check scalars
  qassert(Mesh_get_ndm(mesh) == 1, 'Wrong ndm?');
  qassert(Mesh_get_ndf(mesh) == 1, 'Wrong ndf?');

  % Build the system matrices
  [M,K] = Mesh_assemble_mk(mesh);

  % Check that they look symmetric
  qassert(norm(M-M.','fro')/norm(M,'fro') < 1e-8, 'Unsymm mass');
  qassert(norm(K-K.','fro')/norm(K,'fro') < 1e-8, 'Unsymm stiff');

  % Solve a simple linear field problem
  z = zeros(Mesh_get_numid(mesh),1);
  en = z; en(end) = 1;
  u  = en;
  F  = -K*en;
  u(2:end-1) = K(2:end-1,2:end-1) \ F(2:end-1);
  Mesh_set_u(mesh,u);
  for k = 1:Mesh_numnp(mesh)
    if id(k)
      qassert(abs(Mesh_x(mesh,1,k)/3-u(id(k))) < 1e-8, ...
              '1D linear field is wrong!');
    end
  end

  % Test to see if residual is correct
  R = Mesh_assemble_R(mesh);
  qassert(norm(K*u-R) < 1e-12 * norm(R), 'Mismatch in static residual');

  % Make harmonic and test residual again
  omegat = 1+rand(1);
  Mesh_make_harmonic(mesh, omegat);
  R = Mesh_assemble_R(mesh);
  qassert(norm(K*u - omegat^2*M*u - R) < 1e-12 * norm(R), ...
          'Mismatch in dynamic residual');

  % Clean up
  Mesh_delete(mesh);

end
  

