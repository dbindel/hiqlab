#include <math.h>
#include <algorithm>

#include "material_model.h"


void schur_complement_2D_stress(double *Dpp, double *D);
void schur_complement_2D_stress_te(double *Db, double *D, double alphat);
void schur_complement_2D_stress_pz(double *Cb, double *Cbfull3D);
void schur_complement_2HD_stress_pz(double *Cb, double *Cbfull3D);


/* Function to calculate the elasticity tensor
   for cubic materials.

   lambda    : Material stiffness parameter
   mu        : Material stiffness parameter
   alpha     : Material stiffness parameter
   a         : Directional vector of axis
   b         : Directional vector of axis
   c         : Directional vector of axis
               calculater for a cross b

   v         : Thermal stiffness for 3D
                           Components = {11,22,33,12,23,31}
   vp        :         "         for plane strain(P_components)
   vo        :         "         for plane strain(O_components)
   Vp        :         "         for plane stress(P_components)
   va        :         "         for axis-symm   (A_components)
                Modified from vp
                           Components = {11,22,12}
   C         : Elasticity tensor for 3D
   Ce        :         "         for plane strain(EPSILON)
                           Components = {11,22,12}
   Cs        :         "         for plane stress(SIGMA)
                           Components = {11,22,12}
   Ca        :         "         for axis-symmetric(AXIS)
                           Components = {11,22,33,12}
   Coo       :         "         component for 'o' dofs
   Cooi      : Inverse of Coo
   voTCooivo : vo'*Cooi*vo
   Cb        : Bordered stiffness for Thermoelasticity
               (*:_,e,s,a)
               [ sigma_*    ] = [ C*   -v*] [epsilon_*]
               [beta*tr(eps)]   [ v*'   d ] [theta    ]

   Remarks: p_index is for the PLANE      components {11,22,12}
            o_index is for the OUTOFPLANE components {33,23,31}
            a_index is for the AXISSYM    components {11,22,33,12}
            Otherwise if not noted        components {11,22,33,12,23,31}
*/

void voigt_add_c_aaaa(double c, double* aaaa, double* a);
void voigt_add_c_aaaa(double c, double* aaaa, double* a1,
                      double* a2, double* a3, double* a4);
void voigt_add_c_aaa(double c, double* aaa, double* a1, double* a2, double* a3);
int c_lu_schur(double* A, int ldA, int* ipiv, int n1, int m, int n);


void elasticity_3D(double *C, double lambda, double mu)
{
    std::fill(C, C+6*6, 0);

    // Add 2 * mu * I_rank4
    // Using 2sigma = mu * 2epsilon
    for (int i = 0; i < 3; i++) {
        C[6 *  i    + i  ] += 2*mu;
        C[6 * (i+3) + i+3] +=   mu;
    }

    // Add lambda * I_r2 * I_r2
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++)
            C[6 * i + j] += lambda;
}


void cubic_elasticity_3D(double *C, double lambda, double mu,
                         double alpha, double *a, double *b)
{
    elasticity_3D(C, lambda, mu);

    double c[3] = {a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0]};

    // Add (alpha-2 mu-lambda) * [a*a*a*a + b*b*b*b + c*c*c*c]
    double cc = alpha - lambda - 2 * mu;
    voigt_add_c_aaaa(cc, C, a);
    voigt_add_c_aaaa(cc, C, b);
    voigt_add_c_aaaa(cc, C, c);

}


void hex_elasticity_3D(double *C, double *coeff, double *a, double *b)
{
    double lambda, mu, cc[3];
    lambda = coeff[1];
    mu     =(coeff[0] - coeff[1])/2.0;
    cc[0]  = coeff[3] - coeff[0];
    cc[1]  = coeff[2] - coeff[1];
    cc[2]  = coeff[4] - mu;

    elasticity_3D(C, lambda, mu);

    double c[3] = {a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0]};

    // Add  c[0] * [c*c*c*c]
    //     +c[1] * [a*a*c*c + c*c*a*a + b*b*c*c + c*c*b*b]
    //     +c[2] * [b*c*b*c + b*c*c*b + c*b*b*c + c*b*c*b
    //              c*a*c*a + c*a*a*c + a*c*c*a + a*c*a*c]
    voigt_add_c_aaaa(cc[0], C, c, c, c, c);
    voigt_add_c_aaaa(cc[1], C, a, a, c, c);
    voigt_add_c_aaaa(cc[1], C, c, c, a, a);
    voigt_add_c_aaaa(cc[1], C, b, b, c, c);
    voigt_add_c_aaaa(cc[1], C, c, c, b, b);
    voigt_add_c_aaaa(cc[2], C, b, c, b, c);
    voigt_add_c_aaaa(cc[2], C, b, c, c, b);
    voigt_add_c_aaaa(cc[2], C, c, b, b, c);
    voigt_add_c_aaaa(cc[2], C, c, b, c, b);
    voigt_add_c_aaaa(cc[2], C, c, a, c, a);
    voigt_add_c_aaaa(cc[2], C, c, a, a, c);
    voigt_add_c_aaaa(cc[2], C, a, c, c, a);
    voigt_add_c_aaaa(cc[2], C, a, c, a, c);

}


void elasticity_2D_strain(double *Ce, double lambda, double mu)
{
    double full3D[36];
    elasticity_3D(full3D, lambda, mu);

    int ind[3] = {0, 1, 3};
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            Ce[3*j+i] = full3D[6*ind[j]+ind[i]];
}


void elasticity_2D_stress(double *Cs, double lambda, double mu)
{
    double full3D[36];
    elasticity_3D(full3D, lambda, mu);
    schur_complement_2D_stress(Cs, full3D);
}


void elasticity_axis(double *Ca, double lambda, double mu)
{
    double full3D[36];
    elasticity_3D(full3D, lambda, mu);
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            Ca[4*j+i] = full3D[6*j+i];
}


void cubic_elasticity_2D_strain(double *Ce, double lambda, double mu,
                                double alpha, double *a, double *b)
{
    double full3D[36];
    cubic_elasticity_3D(full3D, lambda, mu, alpha, a, b);

    int ind[3] = {0, 1, 3};
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            Ce[3*j+i] = full3D[6*ind[j]+ind[i]];
}


void cubic_elasticity_2D_stress(double *Cs, double lambda, double mu,
                                double alpha, double *a, double *b)
{
    double full3D[36];
    cubic_elasticity_3D(full3D, lambda, mu, alpha, a, b);
    schur_complement_2D_stress(Cs, full3D);
}


void hex_elasticity_2D_strain(double *Ce, double *coeff, double *a, double *b)
{
    double full3D[36];
    hex_elasticity_3D(full3D, coeff, a, b);

    int ind[3] = {0, 1, 3};
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            Ce[3*j+i] = full3D[6*ind[j]+ind[i]];
}


void hex_elasticity_2D_stress(double *Cs, double *coeff, double *a, double *b)
{
    double full3D[36];
    hex_elasticity_3D(full3D, coeff, a, b);
    schur_complement_2D_stress(Cs, full3D);
}


/* THERMOELASTICITY MODELS */
void thermoelasticity_3D(double *Cb, double lambda, double mu, double alphat)
{
    double C[36],v[6];
    std::fill(Cb, Cb+7*7, 0);
    elasticity_3D(C, lambda, mu);
    v[0] = (3 * lambda + 2 * mu) * alphat;
    v[1] = (3 * lambda + 2 * mu) * alphat;
    v[2] = (3 * lambda + 2 * mu) * alphat;
    v[3] = 0;
    v[4] = 0;
    v[5] = 0;

    // Assemble bordered matrix Cb
    //  [ sigma      ] = [ C   -v] [epsilon]
    //  [beta*tr(eps)]   [ v'   0] [theta  ]
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            Cb[i+7*j] = C[i+6*j];
    for(int i = 0; i < 6; ++i) {
        Cb[6+7*i] = -v[i];
        Cb[i+7*6] =  v[i];
    }
}


void cubic_thermoelasticity_3D(double *Cb, double lambda, double mu,
                               double alpha, double *a, double *b,
                               double alphat)
{
    double C[36],v[6];
    std::fill(Cb, Cb+7*7, 0);
    cubic_elasticity_3D(C, lambda, mu, alpha, a, b);
    v[0] = (2 * lambda + alpha) * alphat;
    v[1] = (2 * lambda + alpha) * alphat;
    v[2] = (2 * lambda + alpha) * alphat;
    v[3] = 0;
    v[4] = 0;
    v[5] = 0;

    // Assemble bordered matrix Cb
    //  [ sigma      ] = [ C   -v] [epsilon]
    //  [beta*tr(eps)]   [ v'   0] [theta  ]
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 6; ++j)
            Cb[i+7*j] = C[i+6*j];
    for(int i = 0; i < 6; ++i) {
        Cb[6+7*i] = -v[i];
        Cb[i+7*6] =  v[i];
    }

}


void schur_complement_2D_stress_te(double *Cb, double *C, double alphat)
{
    // Given
    //
    //  [ sigma_o     ] = [Coo   Cop  -vo] [epsilon_o]
    //  [ sigma_p     ]   [Cop'  Cpp  -vp] [epsilon_p]
    //  [ beta*tr(eps)]   [ vo'   vp'   0] [theta    ]
    //
    // and assuming sigma_o = 0 (plane strain), form the reduced relation
    //
    //  [ sigma_p     ] = [ Cs  -vr] [epsilon_p]
    //  [ beta*tr(eps)]   [ vr'  d ] [theta    ]
    //
    //  Return the bordered matrix
    //
    //  [  Cs   -vr]
    //  [  vr'   d ]
    //
    // Form the thermoelastic stiffness components in the original order
    double ct = alphat * (C[0]+C[1]+C[2]);
    double v[6] = {ct, ct, ct, 0, 0, 0};

    double Creord[49];
    int ind[6] = {2, 4, 5, 0, 1, 3};  // Indices in o-p order

    // Add the border (v) to C and reorder everything
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            Creord[i+7*j] = C[ind[i]+6*ind[j]];
    for (int j = 0; j < 6; ++j){
        Creord[j+7*6] =  v[ind[j]];
        Creord[6+7*j] = -v[ind[j]];
    }
    Creord[6+7*6] = 0.0;

    // Compute the Schur complement in the larger matrix
    int ipiv[3];
    c_lu_schur(Creord, 7, ipiv, 3, 7, 7);

    // Copy the Schur complement into Cb
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i)
            Cb[i+4*j] = Creord[(i+3)+7*(j+3)];
}


void thermoelasticity_2D_strain(double *Cb, double lambda, double mu,
                                double alphat)
{
    double Ce[9],vp[3];
    std::fill(Cb, Cb+4*4, 0);

    elasticity_2D_strain(Ce, lambda, mu);
    vp[0] = alphat * (2 * mu + 3 * lambda);
    vp[1] = alphat * (2 * mu + 3 * lambda);
    vp[2] = 0.0;

    // Assemble bordered matrix Cb
    //  [ sigma_p    ] = [ Ce   -vp] [epsilon]
    //  [beta*tr(eps)]   [ vp'    0] [theta  ]
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            Cb[i+4*j] = Ce[i+3*j];
    for(int i = 0; i < 3; ++i) {
        Cb[3+4*i] = -vp[i];
        Cb[i+4*3] =  vp[i];
    }
}


void thermoelasticity_2D_stress(double *Cb, double lambda,
                                double mu, double alphat)
{
    double full3D[36];
    elasticity_3D(full3D, lambda, mu);
    schur_complement_2D_stress_te(Cb, full3D, alphat);
}


void thermoelasticity_axis(double *Cb, double lambda, double mu, double alphat)
{
    double Ca[16],va[4];
    std::fill(Cb, Cb+5*5, 0);

    elasticity_axis(Ca, lambda, mu);
    va[0] = alphat * (2 * mu + 3 * lambda);
    va[1] = alphat * (2 * mu + 3 * lambda);
    va[2] = alphat * (2 * mu + 3 * lambda);
    va[3] = 0.0;

    // Assemble bordered matrix Cb
    //  [ sigma_a    ] = [ Ca   -va] [epsilon]
    //  [beta*tr(eps)]   [ va'    0] [theta  ]
    for(int i = 0; i < 4; ++i)
        for(int j = 0; j < 4; ++j)
            Cb[i+5*j] = Ca[i+4*j];
    for(int i = 0; i < 4; ++i) {
        Cb[4+5*i] = -va[i];
        Cb[i+5*4] =  va[i];
    }
}


void cubic_thermoelasticity_2D_strain(double *Cb,
                                      double lambda, double mu,
                                      double alpha, double *a,
                                      double *b, double alphat)
{
    double Ce[9],vp[3];
    std::fill(Cb, Cb+4*4, 0);

    cubic_elasticity_2D_strain(Ce, lambda, mu, alpha, a, b);
    vp[0] = alphat * (alpha + 2 * lambda);
    vp[1] = alphat * (alpha + 2 * lambda);
    vp[2] = 0.0;

    // Assemble bordered matrix Cb
    //  [ sigma_p    ] = [ Ce   -vp] [epsilon]
    //  [beta*tr(eps)]   [ vp'    0] [theta  ]
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j)
            Cb[i+4*j] = Ce[i+3*j];
    for(int i = 0; i < 3; ++i) {
        Cb[3+4*i] = -vp[i];
        Cb[i+4*3] =  vp[i];
    }
}


void cubic_thermoelasticity_2D_stress(double *Cb, double lambda,
                                      double mu, double alpha,
                                      double *a, double *b, double alphat)
{
    double full3D[36];
    cubic_elasticity_3D(full3D, lambda, mu, alpha, a, b);
    schur_complement_2D_stress_te(Cb, full3D, alphat);
}


/* PIEZOELECTRIC-ELASTICITY MODELS */
void piezo_elasticity_3D(double *Cb, double lambda, double mu,
                         double *pz, double kds)
{
    double full3D[36],kds_m[9];
    std::fill(kds_m, kds_m+9, 0);
    kds_m[0] = kds_m[4] = kds_m[8] = kds;
    elasticity_3D(full3D, lambda, mu);
    piezo_elasticity_3D(Cb, full3D, pz, kds_m);
}


void piezo_elasticity_3D(double *Cb, double *Cmech, double *pz, double *kds)
{
    std::fill(Cb, Cb+9*9, 0);
    for(int i = 0; i < 6; ++i)
        for(int j = 0; j < 3; ++j) {
            for(int k = 0; k < 6; ++k) {
                Cb[(6+j) + 9*i]-=Cmech[k + i*6] * pz[k + j*6];
            }
            Cb[i + 9*(6+j)] = -Cb[(6+j) + 9*i];
        }

    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) {
            Cb[(6+j) + 9*(6+i)] = kds[j + 3*i];
            for(int k = 0; k < 6; ++k)
                Cb[(6+j) + 9*(6+i)]-=(Cb[k + (6+i)*9]*pz[k + 6*j]);
        }

    for(int j = 0; j < 6; ++j)
        for(int i = 0; i < 6; ++i)
            Cb[i + j*9] = Cmech[i + j*6];

}


void piezo_elasticity_2D_strain(double *Cb, double lambda, double mu,
                                double *pz, double kds)
{
    double Cmech[36],kds_m[9];
    std::fill(kds_m, kds_m+3*3, 0);
    kds_m[0] = kds_m[4] = kds_m[8] = kds;

    elasticity_3D(Cmech, lambda, mu);
    piezo_elasticity_2D_strain(Cb,Cmech,pz,kds_m);
}


void piezo_elasticity_2D_strain(double *Cb, double *Cmech, double *pz, 
                                double *kds)
{
    double Cbfull3D[81];
    std::fill(Cb, Cb+5*5, 0);

    piezo_elasticity_3D(Cbfull3D, Cmech, pz, kds);

    int ind[5] = {0, 1, 3, 6, 7};
    for (int j = 0; j < 5; ++j)
        for (int i = 0; i < 5; ++i)
            Cb[5*j+i] = Cbfull3D[9*ind[j]+ind[i]];

}


void schur_complement_2D_stress_pz(double *Cb, double *Cbfull3D)
{
    // Given
    //
    //  [ sigma_o     ]   [Coo   -eoo  Cop -eop ] [epsilon_o]
    //  [ edisp_o     ] = [eoo'   koo  epo' kop ] [efield _o]
    //  [ sigma_p     ]   [Cop'  -epo  Cpp -epp ] [epsilon_p]
    //  [ edisp_p     ]   [eop'   kop' epp' kpp ] [efield _p]
    //
    //  where,
    //  C*d^T = e = [ epp  epo ]  a 6-by-3 matrix
    //              [ eop  eoo ]
    //
    //
    // and assuming sigma_o = 0 (plane stress), and edisp_o = 0 (plane charge)
    // form the reduced relation
    //
    //  [ sigma_p     ] = [ Cs   -eps  ] [epsilon_p]
    //  [ edisp_p     ]   [ eps' kdps  ] [efield _p]
    //
    //  Return the bordered matrix
    //
    //  [  Cs    -eps ]
    //  [  eps'  kdps ]
    //
    double Creord[81];
    int ind[9] = {2, 4, 5, 8, 0, 1, 3, 6, 7};  // Indices in o-p order

    // Add the border (v) to C and reorder everything
    for (int j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
            Creord[i+9*j] = Cbfull3D[ind[i]+9*ind[j]];

    // Compute the Schur complement in the larger matrix
    int ipiv[4];
    c_lu_schur(Creord, 9, ipiv, 4, 9, 9);

    // Copy the Schur complement into Cb
    for (int j = 0; j < 5; ++j)
        for (int i = 0; i < 5; ++i)
            Cb[i+5*j] = Creord[(i+4)+9*(j+4)];
}


void piezo_elasticity_2D_stress(double *Cb, double lambda, double mu,
                                double *pz, double kds)
{
    double Cmech[36],kds_m[9];
    std::fill(kds_m, kds_m+3*3, 0);
    kds_m[0] = kds_m[4] = kds_m[8] = kds;

    elasticity_3D(Cmech, lambda, mu);
    piezo_elasticity_2D_stress(Cb,Cmech,pz,kds_m);
}


void piezo_elasticity_2D_stress(double *Cb, double *Cmech, double *pz, 
                                double *kds)
{
    double Cbfull3D[81];
    std::fill(Cb, Cb+5*5, 0);

    piezo_elasticity_3D(Cbfull3D, Cmech, pz, kds);
    schur_complement_2D_stress_pz(Cb, Cbfull3D);

}


void schur_complement_2HD_stress_pz(double *Cb, double *Cbfull3D)
{
    // Given
    //
    //  [ sigma_o     ]   [Coo   -eoo  Cop -eop ] [epsilon_o]
    //  [ edisp_o     ] = [eoo'   koo  epo' kop ] [efield _o]
    //  [ sigma_p     ]   [Cop'  -epo  Cpp -epp ] [epsilon_p]
    //  [ edisp_p     ]   [eop'   kop' epp' kpp ] [efield _p]
    //
    //  [ sigma_o     ]   [Coo   Cop  -eoo   -eop ] [epsilon_o]
    //  [ sigma_p     ]   [Cop'  Cpp  -epo   -epp ] [epsilon_p]
    //  [ edisp_o     ] = [eoo'  epo'  koo    kop ] [efield _o]
    //  [ edisp_p     ]   [eop'  epp'  kop'   kpp ] [efield _p]
    //
    //
    //  where,
    //  C*d^T = e = [ epp  epo ]  a 6-by-3 matrix
    //              [ eop  eoo ]
    //
    //
    // and assuming sigma_o = 0 (plane stress), and efield_p = 0 (perpendicular efield)
    // form the reduced relation
    //
    //  [ sigma_p     ] = [ Cs   -eps  ] [epsilon_p]
    //  [ edisp_o     ]   [ eps' kdps  ] [efield _o]
    //
    //  Return the bordered matrix
    //
    //  [  Cs    -eps ]
    //  [  eps'  kdps ]
    //
    double Creord[81];
    int ind[9] = {2, 4, 5, 0, 1, 3, 8, 6, 7};  // Indices in o-p-o-p order

    // Add the border (v) to C and reorder everything
    for (int j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
            Creord[i+9*j] = Cbfull3D[ind[i]+9*ind[j]];

    // Compute the Schur complement in the larger matrix
    int ipiv[3];
    c_lu_schur(Creord, 9, ipiv, 3, 9, 9);

    // Copy the Schur complement into Cb
    for (int j = 0; j < 4; ++j)
        for (int i = 0; i < 4; ++i) {
            Cb[i+4*j] = Creord[(i+3)+9*(j+3)];
        }
}


void piezo_elasticity_2HD_stress(double *Cb, double lambda, double mu,
                                 double *pz, double kds)
{
    double Cmech[36],kds_m[9];
    std::fill(kds_m, kds_m+3*3, 0);
    kds_m[0] = kds_m[4] = kds_m[8] = kds;

    elasticity_3D(Cmech, lambda, mu);
    piezo_elasticity_2HD_stress(Cb,Cmech,pz,kds_m);
}


void piezo_elasticity_2HD_stress(double *Cb, double *Cmech, double *pz, 
                                 double *kds)
{
    double Cbfull3D[81];
    std::fill(Cb, Cb+4*4, 0);

    piezo_elasticity_3D(Cbfull3D, Cmech, pz, kds);
    schur_complement_2HD_stress_pz(Cb, Cbfull3D);

}


void hex_piezo_3D(double *C, double *coeff, double *a, double *b)
{

    std::fill(C, C+3*6, 0);
    double d16,d33,d31;
    d16 = coeff[0];
    d31 = coeff[1];
    d33 = coeff[2];

    double c[3] = {a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0]};
    // Add  d16 * [ a*c*a + a*a*c + b*c*b + b*b*c ]
    //     +d31 * [ c*a*a + c*b*b ]
    //     +d33 * [ c*c*c ]
    voigt_add_c_aaa( d16, C, a, c, a);
    voigt_add_c_aaa( d16, C, a, a, c);
    voigt_add_c_aaa( d16, C, b, c, b);
    voigt_add_c_aaa( d16, C, b, b, c);
    voigt_add_c_aaa( d31, C, c, a, a);
    voigt_add_c_aaa( d31, C, c, b, b);
    voigt_add_c_aaa( d33, C, c, c, c);
}


void hex_dielectric_3D(double *C, double *coeff, double *a, double *b)
{
    std::fill(C, C+3*3, 0);
    double kds1, kds3;
    kds1 = coeff[0];
    kds3 = coeff[1];

    C[0] = C[4] = C[8] = kds1;

    double c[3] = {a[1]*b[2] - a[2]*b[1],
                   a[2]*b[0] - a[0]*b[2],
                   a[0]*b[1] - a[1]*b[0]};

    // Add (kds3-kds1) * [ c*c ]
    for(int i = 0; i < 3; ++i)
        for(int j = 0; j < 3; ++j) {
            C[3*i+j] += (kds3-kds1) * c[i] * c[j];
        }

}


/* HELPER FUNCTIONS */
void voigt_add_c_aaaa(double c, double* aaaa, double* a)
{
    // Form rank 4 outer product (flattened)
    int ind[6] = {0,1,2,0,1,2};
    int jnd[6] = {0,1,2,1,2,0};

    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            aaaa[6*j+i] += c * a[ind[j]] * a[jnd[j]] * a[ind[i]] * a[jnd[i]];
}


void voigt_add_c_aaaa(double c, double* aaaa, double* a1, double* a2, 
                      double* a3, double* a4)
{
    // Form rank 4 outer product (flattened)
    int ind[6] = {0,1,2,0,1,2};
    int jnd[6] = {0,1,2,1,2,0};

    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            aaaa[6*j+i] += c * a1[ind[j]] * a2[jnd[j]] * a3[ind[i]] * a4[jnd[i]];
}


void voigt_add_c_aaa(double c, double* aaa, double* a1, double* a2, double* a3)
{

    // Form rank 4 outer product (flattened)
    int ind[6] = {0,1,2,0,1,2};
    int jnd[6] = {0,1,2,1,2,0};

    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 6; ++i)
            aaa[6*j+i] += c * a1[ind[j]] * a2[ind[i]] * a3[jnd[i]];
}


void schur_complement_2D_stress(double *Cs, double *C)
{
    // Given
    //
    //  [sigma_o] = [Coo  Cop] [epsilon_o]
    //  [sigma_p]   [Cop' Cpp] [epsilon_p]
    //
    // and assuming sigma_o = 0 (plane stress), form the reduced relation
    //
    //   sigma_p = Cs*epsilon_p

    // Reorder the matrix
    double Creord[36];
    int ind[6] = {2, 4, 5, 0, 1, 3};  // Indices in o-p order
    for (int j = 0; j < 6; ++j)
        for (int i = 0; i < 6; ++i)
            Creord[i+6*j] = C[ind[i]+6*ind[j]];

    // Compute the Schur complement in the larger matrix
    int ipiv[3];
    c_lu_schur(Creord, 6, ipiv, 3, 6, 6);

    // Copy the Schur complement into Cs
    for (int j = 0; j < 3; ++j)
        for (int i = 0; i < 3; ++i)
            Cs[i+3*j] = Creord[(i+3)+6*(j+3)];
}


template<class T>
void swap(T& a, T& b)
{
    T tmp = a;
    a = b;
    b = tmp;
}


int c_lu_schur(double* A, int ldA, int* ipiv, int n1, int m, int n)
{
    for (int j = 0; j < n1; ++j) {

        // Find the pivot and test for singularity
        int jp = j;
        for (int i = j+1; i < n1; ++i)
            if (fabs(A[jp+ldA*j]) < fabs(A[i+ldA*j]))
                jp = i;
        ipiv[j] = jp;
        if (A[ldA*j+jp] == 0)
            return j;

        // Apply the interchange to columns 1:n
        if (jp != j)
            for (int i = 0; i < n; ++i)
                swap(A[j+i*ldA], A[jp+i*ldA]);

        // Compute elements j+1:m of jth column
        for (int i = j+1; i < m; ++i)
            A[i+j*ldA] /= A[j+ldA*j];

        // Update trailing submatrix
        for (int jj = j+1; jj < n; ++jj)
            for (int ii = j+1; ii < m; ++ii)
                A[ldA*jj+ii] -= A[ldA*j+ii] * A[ldA*jj+j];

    }
    return 0;
}
