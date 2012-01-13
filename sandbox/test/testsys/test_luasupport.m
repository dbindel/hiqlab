L = Lua_open;

% -- Test valid sets/gets
Lua_set_string(L, 'stringv', 'a stringv');
Lua_set_double(L, 'doublev', 1.234);
Lua_set_array(L, 'arrayv', [1 2 3 4]);
Lua_dofile(L, 'test_luasupport.lua');
s = Lua_get_string(L, 'stringv');
d = Lua_get_double(L, 'doublev');
a = Lua_get_array(L, 'arrayv');
qassert( strcmp(s, 'a string'), 'Incorrect string value' );
qassert( d == 1.357,            'Incorrect double value' );
qassert( all(a == [1 2 3 5]),   'Incorrect array value');

% -- Test gets of non-existent variables
qassert(isempty(Lua_get_string(L, 'foo')), 'Incorrect string (2)');
% -- This check should be reconsidered since calling function on
%    non-existent double gives error
%qassert(Lua_get_double(L, 'foo') == 0,     'Incorrect double (2)');
qassert(isempty(Lua_get_array(L, 'foo')),  'Incorrect array (2)');

% -- Test string gets of non-strings
qassert(str2num(Lua_get_string(L, 'doublev')) == 1.357, ...
        'Bad number to string');
qassert(isempty(Lua_get_string(L, 'arrayv')), 'Bad array to string');

% -- Test numeric gets of non-numbers
% -- This check should be reconsidered since calling function on
%    non-existent double gives error
%qassert(Lua_get_double(L, 'stringv') == 0, 'Bad string to number');
%qassert(Lua_get_double(L, 'arrayv')  == 0, 'Bad array to number');

% -- Test array gets of non-arrays
qassert(isempty(Lua_get_array(L, 'stringv')), 'Bad string to array');
qassert(isempty(Lua_get_array(L, 'doublev')), 'Bad double to array');

Lua_close(L);
