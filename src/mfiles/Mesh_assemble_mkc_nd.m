function [M,K,C, su,sf] = Mesh_assemble_mkc_nd(mesh)

[M,K,C] = Mesh_assemble_mkc(mesh);
[su,sf, M,K,C] = Mesh_matscale(mesh, M,K,C);
