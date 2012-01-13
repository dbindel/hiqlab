#ifndef ADJSTRUCTURE_H
#define ADJSTRUCTURE_H

#include "mesh.h"

#include <vector>

void build_adj(Mesh* mesh, std::vector<int>& jc, std::vector<int>& ir);
void build_adj(Mesh* mesh, int* jc, int* ir);
int  build_adj_nnz(Mesh* mesh);

#endif /* ADJSTRUCTURE_H */
