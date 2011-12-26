-- Loading the Lua input file
runfile = loadfile('test_clean_mesh.lua')

-- Specify parameters
--order = 1;
--dense = 0.25e-6;

-- Run the Lua input file
runfile()

-- Initialize the mesh
mesh:initialize()
mesh:apply_bc()
numid = mesh:get_numid()
print('--- Untied Mesh info ---')
print('Numid :',numid)
print('Numelt:',mesh:numelt())
print('Numnp :',mesh:numnp())
print('Nen   :',mesh:get_nen())
print('Numdof:',mesh:get_ndf()*mesh:numnp())
print('X     :')
for i = 0,mesh:numnp()-1 do
    io.write('   Node[',i,']:')
    for j = 0,mesh:get_ndm()-1 do
        io.write(mesh:x(j,i),'  ')
    end
    io.write('\n')
end
print('ID    :')
for i = 0,mesh:numnp()-1 do
    io.write('   Node[',i,']:')
    for j = 0,mesh:get_ndf()-1 do
        io.write(mesh:id(j,i),'  ')
    end
    io.write('\n')
end
print('IX    :')
for i = 0,mesh:numelt()-1 do
    io.write('   Elem[',i,']:')
    for j = 0,mesh:get_nen(i)-1 do
        io.write(mesh:ix(j,i),'  ')
    end
    io.write('\n')
end

print('\n')
mesh2 = mesh:get_clean_mesh()
mesh2:set_bc(bc_function)
mesh2:initialize()
mesh2:apply_bc()
print('--- Tied Mesh info ---')
print('Numid :',mesh2:get_numid())
print('Numelt:',mesh2:numelt())
print('Numnp :',mesh2:numnp())
print('Nen   :',mesh2:get_nen())
print('Numdof:',mesh2:get_ndf()*mesh2:numnp())
print('X     :')
for i = 0,mesh2:numnp()-1 do
    io.write('   Node[',i,']:')
    for j = 0,mesh2:get_ndm()-1 do
        io.write(mesh2:x(j,i),'  ')
    end
    io.write('\n')
end
print('ID    :')
for i = 0,mesh2:numnp()-1 do
    io.write('   Node[',i,']:')
    for j = 0,mesh2:get_ndf()-1 do
        io.write(mesh2:id(j,i),'  ')
    end
    io.write('\n')
end
print('IX    :')
for i = 0,mesh2:numelt()-1 do
    io.write('   Elem[',i,']:')
    for j = 0,mesh2:get_nen(i)-1 do
        io.write(mesh2:ix(j,i),'  ')
    end
    io.write('\n')
end

-- Assemble stiffness
K = mesh:assemble_dR(1,0,0)
K2= mesh2:assemble_dR(1,0,0)

-- Assemble the loading vector
QArray_type = 0
F = QArray:new(numid, 1, QArray_type)
F2= QArray:new(numid, 1, QArray_type)
mesh:assemble_R()
mesh2:assemble_R()
mesh:get_f(F)
mesh2:get_f(F2)
F:mul(-1)
F2:mul(-1)

-- Assemble LHS
U = QArray:new(numid, 1, QArray_type)
U2= QArray:new(numid, 1, QArray_type)

-- Set up the linear system and solve
K:solve(U, F)
K2:solve(U2, F2)
mesh:set_u(U)
mesh2:set_u(U2)

-- Check residual
R = QArray:new(numid, 1, QArray_type)
R2= QArray:new(numid, 1, QArray_type)
K:apply(U, R)
K2:apply(U2, R2)
R:sub(F)
R2:sub(F2)
resid = R:normf()
resid2= R2:normf()
print('Norm of residual:',resid)
print('Norm of residual:',resid2)

-- Write file for plotting
---[[
mesh2:set_u(U)
dxf = DXFile:new('cleanmesh2d')
dxf:writemesh(mesh2)
dxf:delete()
--]]

-- Print results
U:print()
U2:print()
F:print()
F2:print()

-- Delete objects
R:delete()
U:delete()
F:delete()
K:delete()
R2:delete()
U2:delete()
F2:delete()
K2:delete()

mesh2:delete()
mesh:delete()
