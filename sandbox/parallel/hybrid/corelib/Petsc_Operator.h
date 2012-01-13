#ifndef _PETSC_OPERATOR_H_
#define _PETSC_OPERATOR_H_

#include <iostream>

#include "Epetra_Map.h"
#include "AztecOO.h"
#include "Epetra_LinearProblem.h"

#include "trilinos_epetra_vector.h"
#include "trilinos_epetra_matrix.h"
#include "trilinos_epetra_operator.h"
#include "trilinos_epetra_linearproblem.h"
#include "trilinos_super_matrix.h"
#include "trilinos_super_vector.h"
#include "Epetra_MpiComm.h"

#include "petscksp.h"

#include "mesh.h"

class Epetra_BlockMap;
class Epetra_Comm;


/** Petsc_Operator: An implementation of the 
 *  Petsc class. 
 */

class Petsc_Operator: public virtual Epetra_Operator {
      
 public:

    /** Constructor. */
    Petsc_Operator::Petsc_Operator(KSP ksp, Epetra_Map* dm, Epetra_Map* rm,
		    Epetra_CrsMatrix* B);
    Petsc_Operator::Petsc_Operator(KSP ksp, Epetra_Map* dm,
		    Epetra_CrsMatrix* B);
    Petsc_Operator(KSP ksp, Epetra_CrsMatrix* B);
    virtual ~Petsc_Operator();
  
    /** Multiply by the transpose  
     * (This method has no effect and returns -1 as error code.) 
     */
    int SetUseTranspose(bool UseTranspose){return(-1);};

    /** Apply inverse to a Epetra_MultiVector_Complex X
     *  (This method has no effect and returns -1 as error code.)
     */
    int ApplyInverse(const Epetra_MultiVector& X, 
		           Epetra_MultiVector& Y) const{return 0;};
    
    /** Apply to a Epetra_MultiVector_Complex X
     */
    int Apply(const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const;

    /** Returns the infinity norm of the global matrix. 
     *warning This method must not be called unless HasNormInf() returns true.
     */
    double NormInf() const {return(0.0);};

    /** Returns a character string describing the operator */
    const char * Label() const {return(Label_);};

    /** Returns the current UseTranspose setting. */
    bool UseTranspose() const {return(false);};
  
    /** Returns true if this object can provide an approximate Inf-norm, false otherwise. */
    bool HasNormInf() const{return(false);};
  
    /** Returns a pointer to the Epetra_Comm communicator associated 
     * with this operator. 
     */
    const Epetra_Comm & Comm() const{return(dmap->Comm());};
  
    /** Returns the Epetra_BlockMap object associated with the domain of 
     * this matrix operator. 
     */
    const Epetra_Map & OperatorDomainMap() const
    {
      if (UseTranspose()) return(*dmap);
      else return(*rmap);
    }
  
    /** Returns the Epetra_BlockMap object associated with the range of 
     * this matrix operator. 
     */
    const Epetra_Map & OperatorRangeMap() const
    {
      if (UseTranspose()) return(*rmap);
      else return(*dmap);
    }

    /** Sets mesh object which is used for constructing maps or
     *  PCSetCoordinates
     */
    void SetMesh(Mesh* m);
    void SetFromOptions();

 private:
    void SetReducedMap();
    void set_petsc_map();
    void reduced2full(Epetra_MultiVector* xtmp,
                                  Epetra_MultiVector* xtmp_f) const;
    void full2reduced(Epetra_MultiVector* xtmp_f,
                                  Epetra_MultiVector* xtmp) const;

 protected:

    KSP ksp;
    Epetra_CrsMatrix* M;
    Epetra_Map* dmap;
    Epetra_Map* rmap;
    Epetra_Map* petsc_map;
    Epetra_Map* petsc_map_r;
    int is_petsc_same;    // whether petsc_map and petsc_map_r are same
    Mesh* mesh;
    char * Label_;
};

void Trilinos2PetscMultiVector(Epetra_MultiVector* xtmp_f, Vec* rhs);
void Petsc2TrilinosMultiVector(Vec* lhs, Epetra_MultiVector* xtmp_f);

#endif /* _PETSC_OPERATOR_H_ */

