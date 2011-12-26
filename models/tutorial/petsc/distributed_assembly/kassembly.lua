-- Loading the Lua input file
require 'meshes1.lua'

-- Read in a partition
nmesh = MPI_Comm_rank(PETSC_COMM_WORLD);
mesh  = Mesh_Partition:new(nmesh,2)

mpr = Mesh_Partition_Reader:new()
mesh:set_partition(mpr,'mesh',mesh_func)

-- Assemble stiffness
mat_type = 2
is_reduced=1
Kz = Mesh_Partition_assemble_dR_petsc(mesh,1,0,0,mat_type,is_reduced)

-- Assemble rhs
Fz = Mesh_Partition_assemble_R_petsc(mesh,is_reduced)

-- print results
--pv = PETSC_VIEWER_STDOUT_WORLD
pv = PetscViewerASCIIOpen('petscmat.m')
PetscViewerSetFormat(pv,'matlab')
MatView(Kz,pv)
VecView(Fz,pv)
PetscViewerDestroy(pv)

--[[
print('\n\n')
print('--- Mesh entire info ---')
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

-- Clean up reading environment
MatDestroy(Kz)
VecDestroy(Fz)
mpr:delete()
mesh:delete()
