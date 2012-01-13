
[mesh,L] = Mesh_load('beammesh.lua');

% -- Check that M and K are assembled okay and can transfer
[M,K] = Mesh_assemble_mk(mesh);
Lua_dostring(L, 'M1 = mesh:assemble_dR(0,0,1)');
Lua_dostring(L, 'K1 = mesh:assemble_dR(1,0,0)');
M1 = Lua_get_array(L, 'M1');
K1 = Lua_get_array(L, 'K1');
qassert(norm(M1-M,1) == 0, 'M does not match');
qassert(norm(K1-K,1) == 0, 'K does not match');

% -- Compute MF and K\F
F = Mesh_get_sense_f(mesh, 'force_tip');
MF0 = M*F;
KF0 = K\F;

% -- Try to do the same computation with the default QArray types
Lua_set_array(L, 'F1', F);
Lua_dostring(L, 'MF1 = M1:apply(F1)');
Lua_dostring(L, 'KF1 = K1:solve(F1)');
MF1 = Lua_get_array(L, 'MF1');
KF1 = Lua_get_array(L, 'KF1');
qassert(norm(KF1-KF0)/norm(KF0) < 1e-8, 'Mismatch applying K^-1');
qassert(norm(MF1-MF0)/norm(MF0) < 1e-8, 'Mismatch applying M'   );

% -- Try to do the same computation with explicitly allocated QArrays
Lua_dostring(L, 'MF2 = QArray:new(M1:get_n(),1)');
Lua_dostring(L, 'KF2 = QArray:new(K1:get_n(),1)');
Lua_dostring(L, 'M1:apply(F1,MF2)');
Lua_dostring(L, 'K1:solve(KF2,F1)');
MF2 = Lua_get_array(L, 'MF2');
KF2 = Lua_get_array(L, 'KF2');
qassert(norm(KF1-KF2) == 0, 'Mismatch between two solve calls');
qassert(norm(MF1-MF2) == 0, 'Mismatch between two apply calls');

% -- Try to do the same computation with a complex QArray
Lua_dostring(L, 'F3  = QArray:new(F1:m(),F1:n(),1)');
Lua_dostring(L, 'F3:copy(F1)');
Lua_dostring(L, 'MF3 = QArray:new(M1:get_n(),1,1)');
Lua_dostring(L, 'KF3 = QArray:new(K1:get_n(),1,1)');
Lua_dostring(L, 'M1:apply(F3,MF3)');
Lua_dostring(L, 'K1:solve(KF3,F3)');
MF3 = Lua_get_array(L, 'MF3');
KF3 = Lua_get_array(L, 'KF3');
qassert(norm(KF3-KF0)/norm(KF0) < 1e-8, 'Mismatch applying K^-1');
qassert(norm(MF3-MF0)/norm(MF0) < 1e-8, 'Mismatch applying M'   );

% -- Try to do the same computation with the other complex QArray
Lua_dostring(L, 'F4  = QArray:new(F1:m(),F1:n(),2)');
Lua_dostring(L, 'F4:copy(F1)');
Lua_dostring(L, 'MF4 = QArray:new(M1:get_n(),1,2)');
Lua_dostring(L, 'KF4 = QArray:new(K1:get_n(),1,2)');
Lua_dostring(L, 'M1:apply(F4,MF4)');
Lua_dostring(L, 'K1:solve(KF4,F4)');
MF4 = Lua_get_array(L, 'MF4');
KF4 = Lua_get_array(L, 'KF4');
qassert(norm(KF4-KF0)/norm(KF0) < 1e-8, 'Mismatch applying K^-1');

Mesh_delete(mesh);
