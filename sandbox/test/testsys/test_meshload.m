% Test to see if Mesh_load correctly transfers information

p = [];
p.s = 'A string';
p.n = 1.5;
p.a = [4 3 2 1];
p.br = 1:5;
p.bc = (1:5)';

[mesh,L] = Mesh_load('testmesh.lua', p);

s = Lua_get_string(L, 's');
n = Lua_get_double(L, 'n');
a = Lua_get_array(L, 'a');
a2 = Lua_get_array(L, 'a2');
br = Lua_get_array(L, 'br');
bc = Lua_get_array(L, 'bc');

qassert(strcmp(s, 'A STRING'), 'String was not correctly transformed');
qassert(n == 0.5,              'Double was not correctly transformed');
qassert(all(a == [5 4 3 2]),   'Array was not correctly transformed');
qassert(all(a2 == [4 3 2]),    'Array was not correctly subscripted');
qassert(all(br == (1:5).^2),       'Row was not correctly subscripted');
qassert(all(bc == ((1:5)').^2),    'Col was not correctly subscripted');

Mesh_delete(mesh);
