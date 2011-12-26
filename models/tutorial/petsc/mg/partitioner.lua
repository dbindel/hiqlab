-- Loading the Lua input file
require 'meshes1.lua'

-- Write out to partition info
mnpw = Mesh_Partition_Writer:new()

-- Construct fine(to) mesh
dense_x=  dense_xa[nlevels]
dense_y=  dense_ya[nlevels]
mesh_t = Mesh_Add_Block:new(ndm)
mesh_t = mesh_func(mesh_t)
--mesh_t:initialize()
mesh_t:initialize_minimal()

-- Serial partitioner mesh_t
meshp = Mesh_Partitioner_METIS:new(nparts,mesh_t)
meshp:partition()
mnpw:write_mesh(meshp,fnamea[nlevels])

-- Go through levels
for il = nlevels-1,1,-1 do

    -- Construct coarse(from) mesh
    dense_x=  dense_xa[il]
    dense_y=  dense_ya[il]
    mesh_f = Mesh_Add_Block:new(ndm)
    mesh_f = mesh_func(mesh_f)
    --mesh_f:initialize()
    mesh_f:initialize_minimal()

    -- Construct mapping between meshes
    mnpw:write_nemap(meshp,mesh_f,mesh_t,fpnamea[il])

    -- Extract necessary data and delete previous data
    pepair_data = PEPair_Data:new(meshp,mesh_f,mesh_t)
    meshp:delete()
    mesh_t:delete()

    -- Serial partitioner mesh_f conforming to partition mesh_t
    meshp = Mesh_Partitioner_Conform:new(nparts,mesh_f)
    meshp:build_pepairs_data(pepair_data)
    pepair_data:delete()
    meshp:partition()
    mnpw:write_mesh(meshp,fnamea[il])

    -- Shift the FROM mesh to the TO mesh
    mesh_t = mesh_f

end

-- clean up
meshp:delete()
mesh_t:delete()
mnpw:delete()
