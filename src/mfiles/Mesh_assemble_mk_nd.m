function [M,K, su,sf] = Mesh_assemble_mk_nd(mesh)

[M,K] = Mesh_assemble_mk(mesh);
[su,sf, M,K] = Mesh_matscale(mesh, M,K);
