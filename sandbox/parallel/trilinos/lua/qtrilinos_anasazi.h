#ifndef _QTRILINOS_ANASAZI_H
#define _QTRILINOS_ANASAZI_H

#include "trilinos_anasazi.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

// -- Make sure this is consistent with Anasazis definition
enum MsgType {
    Anasazi_Error = 0,                  /*!< Errors [ always printed ] */
    Anasazi_Warning = 0x1,              /*!< Internal warnings */
    Anasazi_IterationDetails = 0x2,     /*!< Approximate eigenvalues, errors */
    Anasazi_OrthoDetails = 0x4,         /*!< Orthogonalization/orthonormalization details */
    Anasazi_FinalSummary = 0x8,         /*!< Final computational summary */
    Anasazi_TimingDetails = 0x10,       /*!< Timing details */
    Anasazi_Debug = 0x20                /*!< Debugging information */
};

enum SolverType {
    Anasazi_BlockKrylovSchur = 0,
    Anasazi_BlockDavidson    = 1,
    Anasazi_LOBPCG           = 2
};

#endif /* _QTRILINOS_ANASAZI_H */
