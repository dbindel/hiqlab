#include <iostream>
#include <complex>
#include "petsc_mesh.h"
#include "qpassembly_petsc.h"

void Mesh_SetU_Petsc_Vec(Mesh* mesh, Vec v, int is_reduced)
{

    int NumGlobalElements;
    int myid;
    Vec pvec;
    VecScatter scatter;

    MPI_Comm_rank(PETSC_COMM_WORLD,&myid);
    VecGetSize(v,&NumGlobalElements);

    // -- Conduct scatter to zero
    // -- Construct a vector which has all nodes on root
    VecScatterCreateToZero(v,&scatter,&pvec);
    VecScatterBegin(scatter, v, pvec, INSERT_VALUES, SCATTER_FORWARD);
    VecScatterEnd  (scatter, v, pvec, INSERT_VALUES, SCATTER_FORWARD);

    // -- Extract view of the vector and copy
    if (myid==0) {
        int ind[NumGlobalElements];
        PetscScalar* u;
        PetscMalloc(NumGlobalElements*sizeof(PetscScalar),&u);
        for (int i = 0; i < NumGlobalElements; ++i)
            ind[i] = i;
        VecGetValues(pvec, NumGlobalElements, ind, u);

        if (is_reduced==QP_ASSEMBLE_REMOVE_DIRICHLET) {
#ifdef PETSC_USE_COMPLEX
            mesh->set_uz(u);
#else
            mesh->set_u(u);
#endif
        } else if (is_reduced==QP_ASSEMBLE_IDENTITY_DIRICHLET) {
            int numid = mesh->get_numid();
            PetscScalar* ur;
            PetscMalloc(NumGlobalElements*sizeof(numid),&ur);
            for (int i = 0; i < NumGlobalElements; ++i)
                if (mesh->id(i)>=0)
                    ur[mesh->id(i)] = u[i];
#ifdef PETSC_USE_COMPLEX
            mesh->set_uz( ur );
#else
            mesh->set_u( ur );
#endif
            PetscFree(ur);
        }
        PetscFree(u);
    }

    // -- Barrier till copy finishes
    MPI_Barrier(PETSC_COMM_WORLD);

    // -- Clean up
    VecScatterDestroy(scatter);
    VecDestroy(pvec);

}
