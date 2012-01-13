
assert(stringv == 'a stringv', 'Incorrect string value');
assert(doublev == 1.234,       'Incorrect double value');
m,n = arrayv:size()
assert(m == 1 and n == 4, 'Incorrect array dimensions');
for i = 1,4 do
  xr,xi = arrayv(i)
  assert(xr == i and not xi, 'Incorrect element value');
end

stringv = 'a string'
doublev = 1.357
arrayv:set(1,4, 5)

