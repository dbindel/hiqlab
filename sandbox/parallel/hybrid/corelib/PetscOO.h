#ifndef _PETSCOO_H_
#define _PETSCOO_H_

#include <iostream>

#include "petscksp.h"
#include "mesh.h"

#include "qpetsc_pc.h"

class PetscOO  {

  public:

   /** Constructor */
   PetscOO();
   virtual ~PetscOO();

   /** Set Operators */
   int SetOperators(Mat A, Mat P, MatStructure flag);

   /** Set PC type */
   int SetPCType(PCType type);
   int SetPCType(PCType type, Mesh* mesh);

   /** Set option */
   int SetOption(const std::string& iname, const std::string& value);

   /** Solve the problem */
   int Solve(Vec b, Vec x);

   /** Get info */
   double GetResidualNorm();
   int    GetIterationNumber();

  private:

    KSP ksp;
    PC  pc;
    void set_from_options();
};


#endif /* _PETSCOO_H_ */
