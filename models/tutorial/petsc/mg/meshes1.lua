require 'mesh1.lua'

-- Specify parameters
is_removed_unused_nodes = 0

ndm      = 2
nparts   = 2

nlevels  = 3
ordera   = {1, 1, 1}
dense_xa = {2.0e0, 1.0e0, 0.5e0}
dense_ya = {1.0e0, 0.5e0, 0.25e0}
fnamea   = {'mesh_a','mesh_b', 'mesh_c'}
fpnamea  = {'Pab','Pbc'}
