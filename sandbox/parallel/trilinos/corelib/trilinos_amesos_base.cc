#include "hiqlab_trilinos_config.h"
#include "Amesos_config.h"
#include "Amesos.h"

#ifdef HAVE_HIQLAB_SUPERLUDIST
  #ifdef HAVE_AMESOS_SUPERLUDIST
    #include "Amesos_Superludist_Complex.h"
  #endif
#endif

#ifdef HAVE_HIQLAB_MUMPS
  #ifdef HAVE_AMESOS_MUMPS
    #include "Amesos_Mumps_Complex.h"
  #endif
#endif

#include "Epetra_Object.h"

#include "trilinos_amesos_base.h"

static bool verbose = false;

Amesos_BaseSolver_Complex* Amesos_Complex::Create(const char* ClassType,
                          const Epetra_LinearProblem_Complex& LinearProblem )
{
    string CT = ClassType;
    return(Create(CT,LinearProblem));
}

Amesos_BaseSolver_Complex* Amesos_Complex::Create(const string CT,
                          const Epetra_LinearProblem_Complex& LinearProblem )
{

#ifdef HAVE_HIQLAB_SUPERLUDIST
  #ifdef HAVE_AMESOS_SUPERLUDIST
    if ((CT == "Amesos_Superludist") || (CT == "Superludist"))
        return new Amesos_Superludist_Complex(LinearProblem);
  #endif
#endif

#ifdef HAVE_HIQLAB_MUMPS
  #ifdef HAVE_AMESOS_MUMPS
    if ((CT == "Amesos_Mumps") || (CT == "Mumps"))
        return new Amesos_Mumps_Complex(LinearProblem);
  #endif
#endif

    if (verbose) cerr << "Unknown class type:" << CT << endl ;
    return(0);
}

bool Amesos_Complex::Query(const char* ClassType)
{
  string CT = ClassType;
  return(Query(CT));
}

bool Amesos_Complex::Query(const string CT)
{

#ifdef HAVE_HIQLAB_SUPERLU_DIST
  #ifdef HAVE_AMESOS_SUPERLUDIST
    if ((CT == "Amesos_Superludist") || (CT == "Superludist"))
        return true;
  #endif
#endif

#ifdef HAVE_HIQLAB_MUMPS
  #ifdef HAVE_AMESOS_MUMPS
    if ((CT == "Amesos_Mumps") || (CT == "Mumps"))
        return true;
  #endif
#endif

    return(false);
}
