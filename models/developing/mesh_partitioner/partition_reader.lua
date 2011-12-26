-- Loading the Lua input file
require 'meshes1.lua'

-- Prepare Partition_Reader
mnpw = Mesh_Partition_Reader:new()

-- Read in a partition
nmesh = 0;
nlev  = 3;

-- Construct fine(to) mesh
dense_x=  dense_xa[nlev]
dense_y=  dense_ya[nlev]
mesh  = Mesh_Partition:new(nmesh,ndm)
mnpw:read_partition(mesh,fnamea[nlev])
mesh = mesh_func(mesh)
mesh:initialize()
mnpw:read_partition_ids(mesh,fnamea[nlev])

---[[
print('\n\n')
print('--- Mesh_t entire info ---')
mesh_view(mesh)

print('\n\n')
print('GID    :')
for j = 0,mesh:numnp()-1 do
    io.write('   Node[',j,']:')
    local ndf = mesh:get_ndf()
    for k = 0,mesh:get_ndf()-1 do
        io.write(mesh:idg(k+ndf*j),'(',mesh:id(k+ndf*j),')  ')
    end
    io.write('\n')
end
--]]

---[[
dxf = DXFile:new('trial')
dxf:writemesh(mesh)
dxf:delete()
--]]


-- Clean up reading environment
mesh:delete()
mnpw:delete()
