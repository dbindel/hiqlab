#ifndef AMESOS_SUPERLUDIST_COMPLEX_H
#define AMESOS_SUPERLUDIST_COMPLEX_H
#define EPETRA_MPI 1

#include "Amesos_ConfigDefs.h"
#include "Amesos_BaseSolver.h"
#include "Amesos_NoCopiable.h"
#include "Amesos_Utils.h"
#include "Amesos_Time.h"
#include "Amesos_Status.h"
#include "Amesos_Control.h"
#include "Epetra_LinearProblem.h"
#include "Teuchos_RefCountPtr.hpp"
#ifdef EPETRA_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_Comm.h"
#endif

#include "dcomplex.h"
#include "superlu_zdefs.h"
#include "supermatrix.h"
//  SuperLU defines Reduce to be a macro in util.h, this conflicts with Reduce() in Epetra_MultiVector.h
#undef Reduce

#include "trilinos_amesos_base.h"
#include "trilinos_epetra_linearproblem.h"

/** Amesos_Superludist_Complex:  An object-oriented wrapper for Superludist.
 */
class Amesos_Superludist_Complex: public Amesos_BaseSolver_Complex,
                          private Amesos_Time,
                          private Amesos_NoCopiable,
                          private Amesos_Utils,
                          private Amesos_Control,
                          private Amesos_Status
{

 public:

    /** Creates an Amesos_Superludist instance, using an
     *   Epetra_LinearProblem_Complex
     *  Note: The operator in LinearProblem must be an
     *        Epetra_RowMatrix.
     */


    Amesos_Superludist_Complex(const Epetra_LinearProblem_Complex& LinearProblem);

    ~Amesos_Superludist_Complex(void);

    /** Mathematical functions. */
    int SymbolicFactorization() ;
    int NumericFactorization() ;
    int Solve();


    /** Amesos_Superludist does not support transpose at this time.
     *   returns 0 if UseTranspose is set to false, else 1 (failure)
     */
    int SetUseTranspose(bool UseTranspose) { return( UseTranspose?1:0 );};

    /** Get underlying LinearProblem */
    const Epetra_LinearProblem_Complex *GetProblem() const {return(Problem_); };

    /** Returns true if SUPERLUDIST can handle this matrix shape
     *  Returns true if the matrix shape is one that SUPERLUDIST can
     *  handle. SUPERLUDIST only works with square matrices.
     */
    bool MatrixShapeOK() const;

    /** Always returns true. */
    bool UseTranspose() const {return(true);};

    int SetParameters( Teuchos::ParameterList &ParameterList ) ;

    /** Print various timing. */
    void PrintTiming() const;

    /** Print various information about the parameters used by Superludist. */
    void PrintStatus() const;

 private:
    inline const Epetra_Comm& Comm() const
    {
        return(GetProblem()->GetOperator()->Comm());
    };

    inline const Epetra_Import& Importer() const
    {
        return(*(Importer_.get()));
    }

    inline const Epetra_Map& UniformMap() const
    {
        return(*(UniformMap_.get()));
    }

    inline const Epetra_RowMatrix& UniformMatrix() const
    {
        return(*(UniformMatrix_.get()));
    }

    inline const Epetra_RowMatrix& UniformMatrix_i() const
    {
        return(*(UniformMatrixi_.get()));
    }

    inline Epetra_CrsMatrix& CrsUniformMatrix()
    {
        return(*(CrsUniformMatrix_.get()));
    }

    inline Epetra_CrsMatrix& CrsUniformMatrix_i()
    {
        return(*(CrsUniformMatrixi_.get()));
    }

    int RedistributeA();

    int ReFactor();
    int Factor();

    const Epetra_LinearProblem_Complex* Problem_;
    Epetra_RowMatrix *RowMatrixA_  ;  // Problem_->GetOperator()
    Epetra_RowMatrix *RowMatrixAi_ ;  // Problem_->GetOperator_i()

    RefCountPtr<Epetra_Map> UniformMap_;
    RefCountPtr<Epetra_CrsMatrix> CrsUniformMatrix_;
    RefCountPtr<Epetra_CrsMatrix> CrsUniformMatrixi_;
    RefCountPtr<Epetra_RowMatrix> UniformMatrix_;
    RefCountPtr<Epetra_RowMatrix> UniformMatrixi_;
    Teuchos::RefCountPtr<Epetra_Import> Importer_;

    /** Allows FactOption to be used on subsequent calls to pdgssvx
     * from NumericFactorization
     */
    bool ReuseSymbolic_;
    fact_t FactOption_;

    /** redistribute the input matrix prior to calling Superludist */
    bool Redistribute_ ;

    /** True if the SuperLU_DIST's grid has been created (and has to be free'd)
     */
    int GridCreated_ ;
    int FactorizationDone_ ;

    /** True if NumericFactorization() has been successfully called. */
    bool FactorizationOK_ ;

    /** Global dimension of the matrix. */
    int NumGlobalRows_;

    /** Ap, Ai, Aval form the compressed row storage used by SuperLU_DIST */
    vector <int> Ap_;
    vector <int> Ai_;
    vector <doublecomplex> Aval_;
    /** Contains the global ID of local columns. */
    int* Global_Columns_;
    /**  Here are the structures used by Superlu */
    SuperMatrix SuperluA_;
    ScalePermstruct_t ScalePermstruct_;
    LUstruct_t LUstruct_;
    SOLVEstruct_t SOLVEstruct_;

    int nprow_;
    int npcol_;
    /** SuperLU_DIST's grid information. */
    gridinfo_t grid_;
    /** Vector of options */
    superlu_options_t options_;

    bool PrintNonzeros_;
    string ColPerm_;
    string RowPerm_;
    int* perm_c_;
    int* perm_r_;
    string IterRefine_;
    bool ReplaceTinyPivot_;
    bool Equil_;

    int NumNumericFact_;
    int NumSolve_;

};
#endif /* AMESOS_SUPERLUDIST_COMPLEX_H */
