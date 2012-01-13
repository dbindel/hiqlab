
-- Check accesses for a single array

A = QArray:new(5,10)
m,n = A:size()
assert(m == 5 and n == 10, 'Dimensions are incorrect')

for j = 1,10 do
  for i = 1,5 do
    A:set(i,j, i+j)
  end
end

for j = 1,10 do
  for i = 1,5 do
    assert(A(i,j) == i + j)
  end
end

function die_test(i,j) print(A(i,j) == i+j) end

invalid_list = {{0,1}, {1,0},  {0,10}, {5,0},
                {6,1}, {1,11}, {6,10}, {5,11}}

for ii,i in invalid_list do
  assert(not pcall(function(i,j) print(A(i,j)) end, i[1],i[2]), 
         'Uncaught out-of-bounds access')
end


-- Check validity of a subarray view

Aview = A:view(2,3, 2,4)
m,n = Aview:size()
assert(m == 2 and n == 3, 'Incorrect dimensions for Aview')
for i = 1,2 do
  for j = 1,3 do
    assert(Aview(i,j) == i+j+2)
  end
end

invalid_list = {{0,1}, {1,0}, {0,3}, {2,0},
                {3,1}, {1,4}, {3,3}, {2,4}}

for ii,i in invalid_list do
  assert(not pcall(function(i,j) print(Aview(i,j)) end, i[1],i[2]), 
         'Uncaught out-of-bounds access')
end


-- Check to see if make_owner functions properly

for i = 1,2 do
  for j = 1,3 do
    Aview:set(i,j,100)
    assert(Aview(i,j) == 100, 'View assignment went to wrong place')
    assert(A(i+1,j+1) == 100, 'View assignment went to wrong place in A')
    A:set(i+1,j+1, (i+1)+(j+1))
  end
end

Aview:make_owner()
m,n = Aview:size()
assert(m == 2 and n == 3, 'Incorrect dimensions for Aview')

for i = 1,2 do
  for j = 1,3 do
    Aview:set(i,j, 100)
    assert(A(i+1,j+1) == (i+1)+(j+1), 'View assignment went bad')
  end
end


-- Copy left half of A over to the right half
A:view(1,5, 6,10):copy(A:view(1,5, 1,5))
for i = 1,5 do
  for j = 1,5 do
    assert(A(i,j  ) == i+j, 'View copy went bad')
    assert(A(i,j+5) == i+j, 'View copy went bad')
  end
end

-- Make a subcopy

A2 = QArray:new(2,2)
A2:copy(A:view(1,2, 1,2))
for j = 1,5 do
  for i = 1,5 do
    assert(A(i,j) == i+j, 'Error making subcopy')
  end
end


