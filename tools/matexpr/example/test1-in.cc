#include <stdio.h>

void add_plain_stress(double* Ke, double* dN, double E, double nu)
{
    double N_x = dN[0];
    double N_y = dN[1];

    /* <generator matexpr>
      // Displacement formulation of 2D elastic element
      input N_x, N_y;
      input E, nu;
      inout Ke(2,2);

      B = [ N_x, 0; 
            0,   N_y; 
            N_y, N_x ];

      D = E/(1-nu*nu) * 
              [ 1, nu, 0;
                nu, 1, 0;
                0,  0, (1-nu)/2 ];

      Ke += B'*D*B;
    */
}

int main()
{
    double Ke[4] = {0, 0, 0, 0};
    double dN[3] = {1, 2};

    add_plain_stress(Ke, dN, 100, 0.1);
    printf("[%g, %g;\n %g, %g]\n", Ke[0], Ke[2], Ke[1], Ke[3]);
}
