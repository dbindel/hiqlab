#include <cstdio>
#include <cmath>
#include "gaussquad.h"
int main()
{
    if (fabs(gauss_point(0,1)-(0.0000000000000000)) > 1e-15)
        printf("Failure on abscissa (1,1)\n");
    if (fabs(gauss_weight(0,1)-(2.0000000000000000)) > 1e-15)
        printf("Failure on abscissa (1,1)\n");
    if (fabs(gauss_point(0,2)-(-0.5773502691896258)) > 1e-15)
        printf("Failure on abscissa (1,2)\n");
    if (fabs(gauss_weight(0,2)-(0.9999999999999998)) > 1e-15)
        printf("Failure on abscissa (1,2)\n");
    if (fabs(gauss_point(1,2)-(0.5773502691896258)) > 1e-15)
        printf("Failure on abscissa (2,2)\n");
    if (fabs(gauss_weight(1,2)-(0.9999999999999998)) > 1e-15)
        printf("Failure on abscissa (2,2)\n");
    if (fabs(gauss_point(0,3)-(-0.7745966692414835)) > 1e-15)
        printf("Failure on abscissa (1,3)\n");
    if (fabs(gauss_weight(0,3)-(0.5555555555555554)) > 1e-15)
        printf("Failure on abscissa (1,3)\n");
    if (fabs(gauss_point(1,3)-(0.0000000000000000)) > 1e-15)
        printf("Failure on abscissa (2,3)\n");
    if (fabs(gauss_weight(1,3)-(0.8888888888888885)) > 1e-15)
        printf("Failure on abscissa (2,3)\n");
    if (fabs(gauss_point(2,3)-(0.7745966692414834)) > 1e-15)
        printf("Failure on abscissa (3,3)\n");
    if (fabs(gauss_weight(2,3)-(0.5555555555555556)) > 1e-15)
        printf("Failure on abscissa (3,3)\n");
    if (fabs(gauss_point(0,4)-(-0.8611363115940526)) > 1e-15)
        printf("Failure on abscissa (1,4)\n");
    if (fabs(gauss_weight(0,4)-(0.3478548451374541)) > 1e-15)
        printf("Failure on abscissa (1,4)\n");
    if (fabs(gauss_point(1,4)-(-0.3399810435848564)) > 1e-15)
        printf("Failure on abscissa (2,4)\n");
    if (fabs(gauss_weight(1,4)-(0.6521451548625463)) > 1e-15)
        printf("Failure on abscissa (2,4)\n");
    if (fabs(gauss_point(2,4)-(0.3399810435848562)) > 1e-15)
        printf("Failure on abscissa (3,4)\n");
    if (fabs(gauss_weight(2,4)-(0.6521451548625460)) > 1e-15)
        printf("Failure on abscissa (3,4)\n");
    if (fabs(gauss_point(3,4)-(0.8611363115940527)) > 1e-15)
        printf("Failure on abscissa (4,4)\n");
    if (fabs(gauss_weight(3,4)-(0.3478548451374540)) > 1e-15)
        printf("Failure on abscissa (4,4)\n");
    if (fabs(gauss_point(0,5)-(-0.9061798459386639)) > 1e-15)
        printf("Failure on abscissa (1,5)\n");
    if (fabs(gauss_weight(0,5)-(0.2369268850561891)) > 1e-15)
        printf("Failure on abscissa (1,5)\n");
    if (fabs(gauss_point(1,5)-(-0.5384693101056831)) > 1e-15)
        printf("Failure on abscissa (2,5)\n");
    if (fabs(gauss_weight(1,5)-(0.4786286704993667)) > 1e-15)
        printf("Failure on abscissa (2,5)\n");
    if (fabs(gauss_point(2,5)-(0.0000000000000000)) > 1e-15)
        printf("Failure on abscissa (3,5)\n");
    if (fabs(gauss_weight(2,5)-(0.5688888888888889)) > 1e-15)
        printf("Failure on abscissa (3,5)\n");
    if (fabs(gauss_point(3,5)-(0.5384693101056830)) > 1e-15)
        printf("Failure on abscissa (4,5)\n");
    if (fabs(gauss_weight(3,5)-(0.4786286704993669)) > 1e-15)
        printf("Failure on abscissa (4,5)\n");
    if (fabs(gauss_point(4,5)-(0.9061798459386639)) > 1e-15)
        printf("Failure on abscissa (5,5)\n");
    if (fabs(gauss_weight(4,5)-(0.2369268850561889)) > 1e-15)
        printf("Failure on abscissa (5,5)\n");
    if (fabs(gauss_point(0,6)-(-0.9324695142031521)) > 1e-15)
        printf("Failure on abscissa (1,6)\n");
    if (fabs(gauss_weight(0,6)-(0.1713244923791706)) > 1e-15)
        printf("Failure on abscissa (1,6)\n");
    if (fabs(gauss_point(1,6)-(-0.6612093864662645)) > 1e-15)
        printf("Failure on abscissa (2,6)\n");
    if (fabs(gauss_weight(1,6)-(0.3607615730481386)) > 1e-15)
        printf("Failure on abscissa (2,6)\n");
    if (fabs(gauss_point(2,6)-(-0.2386191860831970)) > 1e-15)
        printf("Failure on abscissa (3,6)\n");
    if (fabs(gauss_weight(2,6)-(0.4679139345726903)) > 1e-15)
        printf("Failure on abscissa (3,6)\n");
    if (fabs(gauss_point(3,6)-(0.2386191860831970)) > 1e-15)
        printf("Failure on abscissa (4,6)\n");
    if (fabs(gauss_weight(3,6)-(0.4679139345726918)) > 1e-15)
        printf("Failure on abscissa (4,6)\n");
    if (fabs(gauss_point(4,6)-(0.6612093864662645)) > 1e-15)
        printf("Failure on abscissa (5,6)\n");
    if (fabs(gauss_weight(4,6)-(0.3607615730481386)) > 1e-15)
        printf("Failure on abscissa (5,6)\n");
    if (fabs(gauss_point(5,6)-(0.9324695142031523)) > 1e-15)
        printf("Failure on abscissa (6,6)\n");
    if (fabs(gauss_weight(5,6)-(0.1713244923791704)) > 1e-15)
        printf("Failure on abscissa (6,6)\n");
    if (fabs(gauss_point(0,7)-(-0.9491079123427584)) > 1e-15)
        printf("Failure on abscissa (1,7)\n");
    if (fabs(gauss_weight(0,7)-(0.1294849661688696)) > 1e-15)
        printf("Failure on abscissa (1,7)\n");
    if (fabs(gauss_point(1,7)-(-0.7415311855993945)) > 1e-15)
        printf("Failure on abscissa (2,7)\n");
    if (fabs(gauss_weight(1,7)-(0.2797053914892768)) > 1e-15)
        printf("Failure on abscissa (2,7)\n");
    if (fabs(gauss_point(2,7)-(-0.4058451513773971)) > 1e-15)
        printf("Failure on abscissa (3,7)\n");
    if (fabs(gauss_weight(2,7)-(0.3818300505051189)) > 1e-15)
        printf("Failure on abscissa (3,7)\n");
    if (fabs(gauss_point(3,7)-(-0.0000000000000001)) > 1e-15)
        printf("Failure on abscissa (4,7)\n");
    if (fabs(gauss_weight(3,7)-(0.4179591836734686)) > 1e-15)
        printf("Failure on abscissa (4,7)\n");
    if (fabs(gauss_point(4,7)-(0.4058451513773972)) > 1e-15)
        printf("Failure on abscissa (5,7)\n");
    if (fabs(gauss_weight(4,7)-(0.3818300505051190)) > 1e-15)
        printf("Failure on abscissa (5,7)\n");
    if (fabs(gauss_point(5,7)-(0.7415311855993946)) > 1e-15)
        printf("Failure on abscissa (6,7)\n");
    if (fabs(gauss_weight(5,7)-(0.2797053914892765)) > 1e-15)
        printf("Failure on abscissa (6,7)\n");
    if (fabs(gauss_point(6,7)-(0.9491079123427585)) > 1e-15)
        printf("Failure on abscissa (7,7)\n");
    if (fabs(gauss_weight(6,7)-(0.1294849661688698)) > 1e-15)
        printf("Failure on abscissa (7,7)\n");
    if (fabs(gauss_point(0,8)-(-0.9602898564975361)) > 1e-15)
        printf("Failure on abscissa (1,8)\n");
    if (fabs(gauss_weight(0,8)-(0.1012285362903762)) > 1e-15)
        printf("Failure on abscissa (1,8)\n");
    if (fabs(gauss_point(1,8)-(-0.7966664774136267)) > 1e-15)
        printf("Failure on abscissa (2,8)\n");
    if (fabs(gauss_weight(1,8)-(0.2223810344533746)) > 1e-15)
        printf("Failure on abscissa (2,8)\n");
    if (fabs(gauss_point(2,8)-(-0.5255324099163292)) > 1e-15)
        printf("Failure on abscissa (3,8)\n");
    if (fabs(gauss_weight(2,8)-(0.3137066458778875)) > 1e-15)
        printf("Failure on abscissa (3,8)\n");
    if (fabs(gauss_point(3,8)-(-0.1834346424956498)) > 1e-15)
        printf("Failure on abscissa (4,8)\n");
    if (fabs(gauss_weight(3,8)-(0.3626837833783626)) > 1e-15)
        printf("Failure on abscissa (4,8)\n");
    if (fabs(gauss_point(4,8)-(0.1834346424956500)) > 1e-15)
        printf("Failure on abscissa (5,8)\n");
    if (fabs(gauss_weight(4,8)-(0.3626837833783618)) > 1e-15)
        printf("Failure on abscissa (5,8)\n");
    if (fabs(gauss_point(5,8)-(0.5255324099163288)) > 1e-15)
        printf("Failure on abscissa (6,8)\n");
    if (fabs(gauss_weight(5,8)-(0.3137066458778870)) > 1e-15)
        printf("Failure on abscissa (6,8)\n");
    if (fabs(gauss_point(6,8)-(0.7966664774136270)) > 1e-15)
        printf("Failure on abscissa (7,8)\n");
    if (fabs(gauss_weight(6,8)-(0.2223810344533745)) > 1e-15)
        printf("Failure on abscissa (7,8)\n");
    if (fabs(gauss_point(7,8)-(0.9602898564975364)) > 1e-15)
        printf("Failure on abscissa (8,8)\n");
    if (fabs(gauss_weight(7,8)-(0.1012285362903765)) > 1e-15)
        printf("Failure on abscissa (8,8)\n");
    if (fabs(gauss_point(0,9)-(-0.9681602395076259)) > 1e-15)
        printf("Failure on abscissa (1,9)\n");
    if (fabs(gauss_weight(0,9)-(0.0812743883615744)) > 1e-15)
        printf("Failure on abscissa (1,9)\n");
    if (fabs(gauss_point(1,9)-(-0.8360311073266355)) > 1e-15)
        printf("Failure on abscissa (2,9)\n");
    if (fabs(gauss_weight(1,9)-(0.1806481606948576)) > 1e-15)
        printf("Failure on abscissa (2,9)\n");
    if (fabs(gauss_point(2,9)-(-0.6133714327005905)) > 1e-15)
        printf("Failure on abscissa (3,9)\n");
    if (fabs(gauss_weight(2,9)-(0.2606106964029354)) > 1e-15)
        printf("Failure on abscissa (3,9)\n");
    if (fabs(gauss_point(3,9)-(-0.3242534234038089)) > 1e-15)
        printf("Failure on abscissa (4,9)\n");
    if (fabs(gauss_weight(3,9)-(0.3123470770400032)) > 1e-15)
        printf("Failure on abscissa (4,9)\n");
    if (fabs(gauss_point(4,9)-(-0.0000000000000000)) > 1e-15)
        printf("Failure on abscissa (5,9)\n");
    if (fabs(gauss_weight(4,9)-(0.3302393550012596)) > 1e-15)
        printf("Failure on abscissa (5,9)\n");
    if (fabs(gauss_point(5,9)-(0.3242534234038089)) > 1e-15)
        printf("Failure on abscissa (6,9)\n");
    if (fabs(gauss_weight(5,9)-(0.3123470770400021)) > 1e-15)
        printf("Failure on abscissa (6,9)\n");
    if (fabs(gauss_point(6,9)-(0.6133714327005906)) > 1e-15)
        printf("Failure on abscissa (7,9)\n");
    if (fabs(gauss_weight(6,9)-(0.2606106964029352)) > 1e-15)
        printf("Failure on abscissa (7,9)\n");
    if (fabs(gauss_point(7,9)-(0.8360311073266354)) > 1e-15)
        printf("Failure on abscissa (8,9)\n");
    if (fabs(gauss_weight(7,9)-(0.1806481606948572)) > 1e-15)
        printf("Failure on abscissa (8,9)\n");
    if (fabs(gauss_point(8,9)-(0.9681602395076260)) > 1e-15)
        printf("Failure on abscissa (9,9)\n");
    if (fabs(gauss_weight(8,9)-(0.0812743883615743)) > 1e-15)
        printf("Failure on abscissa (9,9)\n");
    if (fabs(gauss_point(0,10)-(-0.9739065285171719)) > 1e-15)
        printf("Failure on abscissa (1,10)\n");
    if (fabs(gauss_weight(0,10)-(0.0666713443086881)) > 1e-15)
        printf("Failure on abscissa (1,10)\n");
    if (fabs(gauss_point(1,10)-(-0.8650633666889842)) > 1e-15)
        printf("Failure on abscissa (2,10)\n");
    if (fabs(gauss_weight(1,10)-(0.1494513491505809)) > 1e-15)
        printf("Failure on abscissa (2,10)\n");
    if (fabs(gauss_point(2,10)-(-0.6794095682990247)) > 1e-15)
        printf("Failure on abscissa (3,10)\n");
    if (fabs(gauss_weight(2,10)-(0.2190863625159818)) > 1e-15)
        printf("Failure on abscissa (3,10)\n");
    if (fabs(gauss_point(3,10)-(-0.4333953941292472)) > 1e-15)
        printf("Failure on abscissa (4,10)\n");
    if (fabs(gauss_weight(3,10)-(0.2692667193099963)) > 1e-15)
        printf("Failure on abscissa (4,10)\n");
    if (fabs(gauss_point(4,10)-(-0.1488743389816309)) > 1e-15)
        printf("Failure on abscissa (5,10)\n");
    if (fabs(gauss_weight(4,10)-(0.2955242247147533)) > 1e-15)
        printf("Failure on abscissa (5,10)\n");
    if (fabs(gauss_point(5,10)-(0.1488743389816311)) > 1e-15)
        printf("Failure on abscissa (6,10)\n");
    if (fabs(gauss_weight(5,10)-(0.2955242247147526)) > 1e-15)
        printf("Failure on abscissa (6,10)\n");
    if (fabs(gauss_point(6,10)-(0.4333953941292472)) > 1e-15)
        printf("Failure on abscissa (7,10)\n");
    if (fabs(gauss_weight(6,10)-(0.2692667193099961)) > 1e-15)
        printf("Failure on abscissa (7,10)\n");
    if (fabs(gauss_point(7,10)-(0.6794095682990245)) > 1e-15)
        printf("Failure on abscissa (8,10)\n");
    if (fabs(gauss_weight(7,10)-(0.2190863625159813)) > 1e-15)
        printf("Failure on abscissa (8,10)\n");
    if (fabs(gauss_point(8,10)-(0.8650633666889844)) > 1e-15)
        printf("Failure on abscissa (9,10)\n");
    if (fabs(gauss_weight(8,10)-(0.1494513491505807)) > 1e-15)
        printf("Failure on abscissa (9,10)\n");
    if (fabs(gauss_point(9,10)-(0.9739065285171716)) > 1e-15)
        printf("Failure on abscissa (10,10)\n");
    if (fabs(gauss_weight(9,10)-(0.0666713443086884)) > 1e-15)
        printf("Failure on abscissa (10,10)\n");
}
