-- Loading the Lua input file
require 'meshes1.lua'

-- Specify parameters
order  =  ordera[3];
dense_x=  dense_xa[3]
dense_y=  dense_ya[3]

-- Run the Lua input file for single Mesh
mesh = Mesh:new(ndm);
mesh = mesh_func(mesh)
mesh:initialize()

-- Run mesh partitioner
meshp= Mesh_Partitioner_METIS:new(nparts,mesh)
meshp:partition()

-- Construct mesh_manager and meshes on each process and add to Mesh_Manager
mm = Mesh_Manager:new()
mm:create_meshes(meshp,mesh_func)

-- assemble stiffness and rhs with mesh_manager
K = mm:assemble_dR(1,0,0)
mm:assemble_R()
F = QArray:new(mm:compute_numid(),1,2)
mm:get_f(F)

-- assemble stiffness and rhs with mesh
K1= mesh:assemble_dR(1,0,0)
mesh:assemble_R()
F1= QArray:new(mesh:get_numid(),1,2)
mesh:get_f(F1)

-- print results
--[[
K:dump('K.txt')
K:dump('K1.txt')
F:dump('F.txt','%16.16f')
F1:dump('F1.tt','%16.16f')
--]]

--[[
print('\n\n')
print('--- Mesh entire info ---')
mesh_view(mesh)

-- Show info
for i = 1,nparts do
    local cmesh = mm:get_mesh(i-1)
    print('\n\n')
    print('--- Mesh',i,' info ---')
    mesh_view(cmesh)

    print('\n\n')
    print('GID    :')
    for j = 0,cmesh:numnp()-1 do
        io.write('   Node[',j,']:')
        local ndf = cmesh:get_ndf()
        for k = 0,cmesh:get_ndf()-1 do
            io.write(cmesh:idg(k+ndf*j),'(',cmesh:id(k+ndf*j),')  ')
        end
        io.write('\n')
    end
end
--]]

---[[
for i = 1,nparts do
    local cmesh = mm:get_mesh(i-1)
    dxf = DXFile:new(table.concat({'mesh',i}))
    dxf:writemesh(cmesh)
    dxf:delete()
end
--]]


-- clean up
F:delete()
K:delete()
F1:delete()
K1:delete()

mm:delete_meshes()
mm:delete()
meshp:delete()
mesh:delete()
