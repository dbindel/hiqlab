TL_DIR      = $(TRILINOS_DIR)
TL_LIB_DIR  = $(TRILINOS_LIB_DIR)
INCLUDE_TL= -I$(TL_DIR)/packages/amesos/src \
            -I$(TL_DIR)/packages/anasazi/src \
            -I$(TL_DIR)/packages/aztecoo/src \
            -I$(TL_DIR)/packages/epetra/src \
            -I$(TL_DIR)/packages/epetraext/src \
            -I$(TL_DIR)/packages/epetraext/src/inout \
            -I$(TL_DIR)/packages/epetraext/src/transform \
            -I$(TL_DIR)/packages/teuchos/src \
            -I$(TL_DIR)/packages/tritutils/src \
            -I$(TL_DIR)/packages/ifpack/src \
            -I$(TL_DIR)/packages/ml/src/Include \
            -I$(TL_DIR)/packages/ml/src/Utils \
            -I$(TL_DIR)/packages/ml/src/Main \
            -I$(TL_DIR)/packages/ml/src/Krylov \
            -I$(TL_DIR)/packages/ml/src/FEGrid \
            -I$(TL_DIR)/packages/ml/src/Coarsen \
            -I$(TL_DIR)/packages/ml/src/Comm \
            -I$(TL_DIR)/packages/ml/src/MLAPI \
            -I$(TL_DIR)/packages/ml/src/MatrixFree \
            -I$(TL_DIR)/packages/ml/src/NonlinML \
            -I$(TL_DIR)/packages/ml/src/Operator \
            -I$(TL_DIR)/packages/ml/src/Smoother \
            -I$(TL_LIB_DIR)/packages/amesos/src \
            -I$(TL_LIB_DIR)/packages/anasazi/src \
            -I$(TL_LIB_DIR)/packages/aztecoo/src \
            -I$(TL_LIB_DIR)/packages/epetra/src \
            -I$(TL_LIB_DIR)/packages/epetraext/src \
            -I$(TL_LIB_DIR)/packages/teuchos/src \
            -I$(TL_LIB_DIR)/packages/triutils/src \
            -I$(TL_LIB_DIR)/packages/ifpack/src \
            -I$(TL_LIB_DIR)/packages/ml/src

LIB_TL = $(TL_LIB_DIR)/packages/amesos/src/libamesos.a \
         $(TL_LIB_DIR)/packages/anasazi/src/libanasazi.a \
         $(TL_LIB_DIR)/packages/aztecoo/src/libaztecoo.a \
         $(TL_LIB_DIR)/packages/ml/src/libml.a \
         $(TL_LIB_DIR)/packages/ifpack/src/libifpack.a \
         $(TL_LIB_DIR)/packages/epetraext/src/libepetraext.a \
         $(TL_LIB_DIR)/packages/epetra/src/libepetra.a \
         $(TL_LIB_DIR)/packages/teuchos/src/libteuchos.a 
