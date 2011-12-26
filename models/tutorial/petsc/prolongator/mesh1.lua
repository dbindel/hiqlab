-- Include function definition file
require 'common.lua'

function mesh_func(mesh)

-- Remove unused nodes ?
is_removed_unused_nodes = is_removed_unused_nodes or 0

-- Define element type
mtype     = {}
mtype.E   =  1000.0
mtype.nu  =    0.25
local etype = mesh:PMLElastic2d_planestrain(mtype);

-- Define beam geometry
local l =  l or 2
local w =  w or 2

-- Element order, approximate size of elements
order  =  order   or 1
dense_x=  dense_x or 2e0
dense_y=  dense_y or 2e0
meshtol= 1e-6

-- Define coordinates, element connectivity, and addition of
-- nodes and elements to the simultaneously, using the 
-- blocks2d command
mesh:blocks2d({0, l},{-w, 0}, etype, order, dense_x, dense_y)
mesh:blocks2d({0, l},{ 0, w}, etype, order, dense_x, dense_y)

-- Set stretch function
etype:set_stretch(function(x) return 0; end)

-- Tie mesh together
mesh:tie()

if (is_removed_unused_nodes==1) then
    mesh:remove_unused_nodes()
end

-- Define boundary conditions
local function bc_function(x,y)
--    if mesheq(x,0) then return 'u',  0; end;
--    if mesheq(x,l) then return 'f ', 1   ; end;
--    if mesheq(y,-w) then return 'uu', 0, 0; end;
--    if mesheq(y, w) then return 'f ', 1   ; end;
end
mesh:set_bc(bc_function)

return mesh

end

function mesh_view(mesh)

    local numid = mesh:get_numid()
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
end

