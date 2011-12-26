%
% Test harness for Lua-MATLAB interaction facilities in luamatlab.mw
%

L = Lua_open();

% -- Check dofile/dostring and transfer of basic types -- 

Lua_set_string(L, 's', 'foo');
Lua_set_double(L, 'pi2', 2*pi);

Lua_dofile(L, 'test_luamatlab.lua');
Lua_dostring(L, 'pi = pi2/2');

s2 = Lua_get_string(L, 's2');
p  = Lua_get_double(L, 'pi');

assert(strcmp(s2, 'foofoo'));
assert(p == pi);

% -- Check real QArray transfers --

Lua_set_array(L, 'a', [1, 3, 5, 7]);
Lua_dostring(L, 'a = a + 1');
a = Lua_get_array(L, 'a');
assert(norm(a - [2, 4, 6, 8]) == 0);

% -- Check real CSCMatrix transfers --

A = sprand(20,20, 0.2);
x = rand(20,1);
Lua_set_array(L, 'A', A);
Lua_set_array(L, 'x', x);
Lua_dostring(L, 'b = A:apply(x)');
Lua_get_array(L, 'b');
assert(norm(Lua_get_array(L, 'A')-A, 'fro') == 0);
assert(norm(Lua_get_array(L, 'b')-A*x, 'fro') < ...
       1e-14 * norm(A, 'fro') * norm(x, 'fro'));
Lua_dostring(L, 'b:delete()');

% -- Check complex QArray transfers --

Lua_set_array(L, 'a', 1i*[1, 3, 5, 7]);
Lua_dostring(L, 'a = a + 1');
a = Lua_get_array(L, 'a');
assert(norm(a - (1+[1i, 3i, 5i, 7i])) == 0);

% -- Check behavior on nonexistent gets

try
  pass = 1; Lua_get_double(L, 'nonesuch'); pass = 0;
  error('Should not get here!')
catch
  assert(pass == 1);
end

try
  pass = 1; Lua_get_string(L, 'nonesuch'); pass = 0;
  error('Should not get here!')
catch
  assert(pass == 1);
end

try
  pass = 1; Lua_get_array(L, 'nonesuch'); pass = 0;
  error('Should not get here!');
catch
  assert(pass == 1);
end

Lua_close(L);
