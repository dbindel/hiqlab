#include <iostream>
#include "qpetsc_pc.h"
#include "mesh.h"

int PCSetCoordinatesFromMesh(PC pc, Mesh* mesh)
{
#ifdef PETSC_HAVE_PROMETHEUS
    Mat A;
    PetscInt Istart, Iend;
    int ndf = mesh->get_ndf();
    int ndm = mesh->get_ndm();
    int numlnp, snode, enode;
    PetscReal* coords;
//    double* coords;

    // Currently, Prometheus doesn't copy correctly if ndf != 3
    if (mesh == NULL || mesh->get_ndf() != 3)
        return 0;

    // -- Get the operator related with the preconditioning
    PCGetOperators(pc, PETSC_NULL, &A, PETSC_NULL);

    // -- Get the starting and ending indices/nodes
    MatGetOwnershipRange(A, &Istart, &Iend);
    snode = Istart/ndf;
    enode = Iend/ndf;
    coords = new PetscReal[Iend-Istart];
//    coords = new double[Iend-Istart];
    numlnp = (Iend-Istart)/ndf;
    memset(coords, 0, numlnp * ndf * sizeof(double));

    for (int i = snode; i < enode; ++i)
        for (int j = 0; j < mesh->get_ndm(); ++j) {
            coords[(i-snode)*ndm+j] = mesh->x(j,i);
        }
    PetscErrorCode ierr = PCSetCoordinates(pc, ndm, coords);

    // -- Clean up
    delete[] coords;
    return ierr;
#else
    std::cout << "Must configure PETSc with PROMETHEUS\n";
#endif
}
