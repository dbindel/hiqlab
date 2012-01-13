#include "Epetra_MpiComm.h"
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "qtrilinos_epetra.h"

extern Epetra_MpiComm* HIQLAB_Comm;

/* ======== ================ *
 * function CrsMatrix2MATLAB *
 * ======== ================ *
 *
 * Print out a CrsMatrix in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements. Output is finished by "End of Matrix Output".
 *
 *
 * Return code:        true if matrix has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to the ditributed CrsMatrix to
 *                     print out
 */

int CrsMatrix2MATLAB(Epetra_CrsMatrix* Ap )

{
  Epetra_CrsMatrix A = *Ap;
  int MyPID = A.Comm().MyPID();
  int NumProc = A.Comm().NumProc();

  // work only on transformed matrices;
  if( A.IndicesAreLocal() == false ) {
    if( MyPID == 0 ) {
      cerr << "ERROR in "<< __FILE__ << ", line " << __LINE__ << endl;
      cerr << "Function CrsMatrix2MATLAB accepts\n";
      cerr << "transformed matrices ONLY. Please call A.FillComplete()\n";
      cerr << "on your matrix A to that purpose.\n";
      cerr << "Now returning...\n";
    }
    return 0;
  }

  int NumMyRows = A.NumMyRows(); // number of rows on this process
  int NumNzRow;   // number of nonzero elements for each row
  int NumEntries; // number of extracted elements for each row
  int NumGlobalRows; // global dimensio of the problem
  int GlobalRow;  // row in global ordering
  int NumGlobalNonzeros; // global number of nonzero elements

  NumGlobalRows = A.NumGlobalRows();
  NumGlobalNonzeros = A.NumGlobalNonzeros();

  // print out on cout if no filename is provided

  int IndexBase = A.IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1;

  // write on file the dimension of the matrix

  if( MyPID==0 ) {
    cout << "A = spalloc(";
    cout << NumGlobalRows << ',' << NumGlobalRows;
    cout << ',' << NumGlobalNonzeros << ");\n";
  }

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {

    if( MyPID == Proc ) {

      cout << "% On proc " << Proc << ": ";
      cout << NumMyRows << " rows and ";
      cout << A.NumMyNonzeros() << " nonzeros\n";

      // cycle over all local rows to find out nonzero elements
      for( int MyRow=0 ; MyRow<NumMyRows ; ++MyRow ) {

        GlobalRow = A.GRID(MyRow);

        NumNzRow = A.NumMyEntries(MyRow);
        double *Values = new double[NumNzRow];
        int *Indices = new int[NumNzRow];

        A.ExtractMyRowCopy(MyRow, NumNzRow,
                           NumEntries, Values, Indices);
        // print out the elements with MATLAB syntax
        for( int j=0 ; j<NumEntries ; ++j ) {
          cout << "A(" << GlobalRow  + IndexBase
               << "," << A.GCID(Indices[j]) + IndexBase
               << ") = " << Values[j] << ";\n";
        }

        delete Values;
        delete Indices;
      }

    }
    A.Comm().Barrier();
    if( MyPID == 0 ) {
      cout << " %End of Matrix Output\n";
    }
  }

  return 1;

}

/* ======== ============= *
 * function Vector2MATLAB *
 * ======== ============= *
 *
 * Print out a Epetra_Vector in a MATLAB format. Each processor prints out
 * its part, starting from proc 0 to proc NumProc-1. The first line of
 * each processor's output states the number of local rows and of
 * local nonzero elements. Output is finished by "End of Vector Output".
 *
 * Return code:        true if vector has been printed out
 * -----------         false otherwise
 *
 * Parameters:
 * ----------
 *
 * - Epetra_CrsMatrix  reference to vector
 * - Epetra_Map        reference to map
 */

int Vector2MATLAB(Epetra_Vector* vp)
{
  Epetra_Vector v = *vp;

  Epetra_BlockMap Map = v.Map();
  int MyPID = Map.Comm().MyPID();
  int NumProc = Map.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();

  // print out on cout if no filename is provided

  int IndexBase = Map.IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1;

  // write on file the dimension of the matrix

  if( MyPID == 0 ) cout << "v = zeros(" << GlobalLength << ")\n";

  // get update list
  int * MyGlobalElements = Map.MyGlobalElements( );

  int Row;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {

    if( MyPID == Proc ) {

      cout << "% On proc " << Proc << ": ";
      cout << MyLength << " rows of ";
      cout << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        cout << "v(" << MyGlobalElements[Row]+IndexBase
             << ") = " << v[Row] << ";\n";
      }

      if( MyPID == NumProc-1  ) {
        cout << "% End of vector\n";
      }

    }

    Map.Comm().Barrier();
  }

  return 1;

} /* Vector2MATLAB */

/* ======== ============= *
 * function MultiVector2MATLAB *
 * ======== ============= *
 */

int MultiVector2MATLAB(Epetra_MultiVector* vp, int index)
{
  Epetra_Vector v = Epetra_Vector(View, *vp, index);

  Epetra_BlockMap Map = v.Map();
  int MyPID = Map.Comm().MyPID();
  int NumProc = Map.Comm().NumProc();
  int MyLength = v.MyLength();
  int GlobalLength = v.GlobalLength();

  // print out on cout if no filename is provided

  int IndexBase = Map.IndexBase(); // MATLAB start from 0
  if( IndexBase == 0 ) IndexBase = 1;

  // write on file the dimension of the matrix

  if( MyPID == 0 ) cout << "v = zeros(" << GlobalLength << ")\n";

  // get update list
  int * MyGlobalElements = Map.MyGlobalElements( );

  int Row;

  for( int Proc=0 ; Proc<NumProc ; ++Proc ) {

    if( MyPID == Proc ) {

      cout << "% On proc " << Proc << ": ";
      cout << MyLength << " rows of ";
      cout << GlobalLength << " elements\n";

      for( Row=0 ; Row<MyLength ; ++Row ) {
        cout << "v(" << MyGlobalElements[Row]+IndexBase
             << ") = " << v[Row] << ";\n";
      }

      if( MyPID == NumProc-1  ) {
        cout << "% End of vector\n";
      }

    }

    Map.Comm().Barrier();
  }

  return 1;

} /* MultiVector2MATLAB */


Epetra_Vector* qVectorCreate(int n)
{
    Epetra_Map Map(n, 0, *HIQLAB_Comm);
    Epetra_Vector* x = new Epetra_Vector(Map);
    return x;
}

Epetra_MultiVector* qMultiVectorCreate(int m, int n)
{
    Epetra_Map Map(m, 0, *HIQLAB_Comm);
    Epetra_MultiVector* x = new Epetra_MultiVector(Map, n);
    return x;
}

Epetra_MultiVector_Complex* qMultiVectorComplexCreate(int m, int n)
{
    Epetra_Map Map(m, 0, *HIQLAB_Comm);
    Epetra_MultiVector_Complex* x = new Epetra_MultiVector_Complex(Map, n, 0);
    return x;
}

Epetra_Vector_Complex* qVectorComplexCreate(int n)
{
    Epetra_Map Map(n, 0, *HIQLAB_Comm);
    Epetra_Vector_Complex* x = new Epetra_Vector_Complex(Map, 0);
    return x;
}

Epetra_CrsMatrix* qRow2CrsMatrix(Epetra_RowMatrix* ERM)
{
     return dynamic_cast<Epetra_CrsMatrix *>(ERM);
}

Epetra_RowMatrix* qCrs2RowMatrix(Epetra_CrsMatrix* ECM)
{
     return ECM;
}

Epetra_Vector* qMultiVector2Vector(Epetra_MultiVector* EMV)
{
     return dynamic_cast<Epetra_Vector *>(EMV);
}

Epetra_MultiVector* qVector2MultiVector(Epetra_Vector* EV)
{
     return EV;
}

Epetra_Vector_Complex*      qMulti2SingleComplexVector(Epetra_MultiVector_Complex* v)
{
    Epetra_Vector_Complex* evc;

    evc = new Epetra_Vector_Complex(Copy, (*v), (*(v->get_Vz())), 0 );

    return evc;
}

Epetra_MultiVector_Complex* qSingle2MultiComplexVector(Epetra_Vector_Complex* v)
{
    Epetra_MultiVector_Complex* emvc;

    emvc = new Epetra_MultiVector_Complex( (*v), (*(v->get_Vz())) );

    return emvc;
}
