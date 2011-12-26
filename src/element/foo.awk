{
  gsub(/set_local_u[ ]*[(]mesh, /, "set_local_vec(mesh, Mesh::VEC_U, ");
  gsub(/set_local_v[ ]*[(]mesh, /, "set_local_vec(mesh, Mesh::VEC_V, ");
  gsub(/set_local_a[ ]*[(]mesh, /, "set_local_vec(mesh, Mesh::VEC_A, ");
  print
}
