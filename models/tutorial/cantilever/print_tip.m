function print_tip(mesh, L, ytip)

% -- Obtain tip displacement(MATLAB)
p         = Mesh_get_x(mesh);
beam_l    = Lua_get_double(L,'l');
tip_node  = find( (p(1,:)==beam_l) & (p(2,:)==ytip) );
dsp       = Mesh_get_disp(mesh);
tip_disp1 = dsp(2,tip_node);

% -- Obtain tip displacement through Mesh_get_vector
U         = Mesh_get_u(mesh);
sense_pat = Mesh_get_vector(mesh,'tip_displacement2');
tip_disp2 = sense_pat'*U;

% -- Print results
disp(' ');
disp('Method 1: General but tedious');
fprintf('Tip displacement y: %e\n',tip_disp1);
disp(' ');
disp('Method 2: Using Mesh_get_vector');
fprintf('Tip displacement y: %e\n',tip_disp2);

