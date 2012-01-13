
s = string.upper(s)
n = n-1
for i=1,a:n() do
  a:set(1,i, a(i)+1)
end
a2 = a:view(1,1,2,4)

for i=1,5 do
  br[i] = br(i)*i
  bc[i] = bc(i)*i
end

mesh = Mesh:new(1);
