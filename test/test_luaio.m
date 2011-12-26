disp('Testing I/O redirection: ignore output between dashes');
disp('----');
L = Lua_open;

% -- Flush file
fid = fopen('test_luaio.out', 'w');
fclose(fid);

diary('test_luaio.out');
Lua_dostring(L, 'print ''Test successful''');
diary off

Lua_close(L);
disp('----');

fid = fopen('test_luaio.out', 'r');
tline1 = fgetl(fid);
tline2 = fgetl(fid);
fclose(fid);

assert(ischar(tline1) & ~ischar(tline2));
assert(strcmp(tline1, 'Test successful'));
