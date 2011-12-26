#include <stdio.h>
#include <math.h>
#include <complex>

int test_status = 0;


void tassert(int condition, const char* msg)
{
    if (!condition) {
        printf("FAILED: %s\n", msg);
	test_status = -1;
    }
}


int approx(double a, double b, double rtol)
{
    return ( fabs(a-b) < rtol * fabs(a) );
}


/*
 * Check:
 *   Assignments
 *   Addto
 *   Declarations (array, scalar, vector)
 *   Scalar interactions, vector interactions
 *   Negation
 */

void test_scalar_mat()
{
    double add_sv[2], add_vs[2];
    double sub_sv[2], sub_vs[2];
    double mul_sv[2], mul_vs[2];
    double div_vs[2];
    double neg_s;
    double neg_v[2];

    /* <generator matexpr>
    // Test scalar op matrix calculations
    output add_sv(2), add_vs(2);
    output sub_sv(2), sub_vs(2);
    output mul_sv(2), mul_vs(2);
    output div_vs(2);
    output neg_s, neg_v(2);

    v = [1.5; 6];
    s = 3;

    add_sv = s + v;
    add_vs = v + s;
    sub_sv = s - v;
    sub_vs = v - s;
    mul_sv = s * v;
    mul_vs = v * s;
    div_vs = v / s;

    neg_s = -s;
    neg_v = -v;
    */

    tassert(add_sv[0] ==  4.5 && add_sv[1] ==  9.0, "scalar + vector");
    tassert(add_vs[0] ==  4.5 && add_vs[1] ==  9.0, "vector + scalar");
    tassert(sub_sv[0] ==  1.5 && sub_sv[1] == -3.0, "scalar - vectpr");
    tassert(sub_vs[0] == -1.5 && sub_vs[1] ==  3.0, "vector - scalar");

    tassert(mul_sv[0] ==  4.5 && mul_sv[1] == 18.0, "scalar * vector");
    tassert(mul_vs[0] ==  4.5 && mul_vs[1] == 18.0, "vector * scalar");
    tassert(div_vs[0] ==  0.5 && div_vs[1] ==  2.0, "vector / scalar");

    tassert(neg_s == -3,                            "-scalar");
    tassert(neg_v[0] == -1.5 && neg_v[1] == -6,     "-vector");
}


void test_mat_mat()
{
    double add_vv[2];
    double sub_vv[2];
    double mul_mv[2];

    /* <generator matexpr>
    // Test scalar op matrix calculations
    output add_vv(2), sub_vv(2), mul_mv(2);

    A = [ 1, -1;
         -2,  2];
    v1 = [1;  3];
    v2 = [5;  8];

    add_vv = v2 + v1;
    sub_vv = v2 - v1;
    mul_mv = A*v1;
    */

    tassert(add_vv[0] == 6 && add_vv[1] == 11, "vector + vector");
    tassert(sub_vv[0] == 4 && sub_vv[1] ==  5, "vector - vector");
    tassert(mul_mv[0] ==-2 && mul_mv[1] ==  4, "matrix * vector");
}


void test_mat()
{
    double v1, v2, v3, v4;
    double v1b, v2b, v3b, v4b;
    double I[4];

    /* <generator matexpr>
    output v1, v2, v3, v4;
    output v1b, v2b, v3b, v4b;
    output I(2,2);

    A = [1, 2; 3, 4];
    e1 = [1; 0];
    e2 = [0; 1];

    v1 = e1'*A*e1;
    v2 = e1'*A*e2;
    v3 = e2'*A*e1;
    v4 = e2'*A*e2;

    v1b = e1'*(A*e1);
    v2b = e1'*(A*e2);
    v3b = e2'*(A*e1);
    v4b = e2'*(A*e2);

    I = [e1, e2];
    */

    tassert(v1 == 1 && v1b == 1 && 
            v2 == 2 && v2b == 2 &&
            v3 == 3 && v3b == 3, "Matrix multiply");
    tassert(I[0] == 1 && I[1] == 0 && I[2] == 0 && I[3] == 1, "Matrix concat");
}


void test_opt()
{
    double a = 1;
    double b = 2;
    double z1, z2, z3, z4;
    double d1, d2, d3, d4, d5, d6, d7, d8;

    /* <generator matexpr>
    input a, b;
    output z1, z2, z3, z4;
    output d1, d2, d3, d4, d5, d6, d7, d8;

    // Zero removal
    z1 = 0 * b;
    z2 = b * 0;
    z3 = 0 / b;
    z4 = -0;

    // Identity simplification
    d1 = a + 0;
    d2 = 0 + a;
    d3 = a - 0;
    d4 = 0 - a;
    d5 = a * 1;
    d6 = 1 * a;
    d7 = a / 1;

    // Copy prop
    d8 = d1 + z1 + z2 + z3 + z4;

    // Dead code elimination
    c = a + b;

    */

    tassert(z1 == 0 && z2 == 0 && z3 == 0 && z4 == 0, "Zero folding opt");
    tassert(d1 == a && d2 == a && d3 == a && d4 == -a &&
            d5 == a && d6 == a && d7 == a && d8 == a, "Identity folding opt");
}


void test_call()
{
    double four = 4;
    double two;
    double I2[4];

    /* <generator matexpr>
    input four;
    output two = sqrt(four);
    output I2(2,2) = eye(2);
    */

    tassert(two == 2, "Function call");
    tassert(I2[0] == 1 && I2[1] == 0 && I2[2] == 0 && I2[3] == 1, "eye call");
}


void test_subscript()
{
    double A[4] = {1, 2, 3, 4};
    double b[2] = {5, 6};
    double A12, b2;

    /* <generator matexpr>
    input A(2,2), b(2);
    output A12, b2;
    A12 = A(1,2);
    b2 = b(2);
    */

   tassert(A12 == 3, "Matrix subscript");
   tassert(b2  == 6, "Vector subscript");
}


void test_inout()
{
    double a = 1;
    double b[2] = {2, 3};

    /* <generator matexpr>
    inout a, b(2);

    a += 10;
    b += 10;
    */

    tassert(a == 11, "Inout / += on scalar");
    tassert(b[0] == 12 && b[1] == 13, "Inout / += on vector");
}


void test_iconst()
{
    double half;
    /* <generator matexpr>
    output half;
    half = 1/2;
    */

    tassert(half == 0.5, "Integer constant management");
}


void test_range()
{

    double range[3];
    /* <generator matexpr>
    output range(3);
    range = 3:5;
    */

    tassert(range[0] == 3 && range[1] == 4 && range[2] == 5,
            "Range generation");
}


void test_deriv()
{
    double x = 0.4;
    double y = 0.7;
    double Df1[4], Df3[2], Df4[2];
    double Df2, Df5, Df6, Df7, Df8, Df9, Df10;
    /* <generator matexpr>
    input x, y;
    output Df1(2,2), Df2, Df3(2), Df4(2), Df5, Df6, Df7, Df8, Df9, Df10;
    f1 = [x; y];
    Df1 = deriv(f1, [x,y]);
    Df2 = deriv(2*x, x);
    Df3 = deriv(f1'*f1, [x;y]);
    Df4 = deriv(x-y, [x;y]);
    Df5 = deriv(-y, y);
    Df6 = deriv(1/x, x);
    Df7 = deriv(sqrt(x*x), x);
    Df8 = deriv(log(x), x);
    Df9 = deriv(log(exp(x)), x);
    Df10 = deriv(exp(x*y), x);
    */
    tassert(Df1[0] == 1 && Df1[1] == 0 && Df1[2] == 0 && Df1[3] == 1,
            "Basic differentiation");
    tassert(Df2 == 2,                       "Differentiate scalar multiply");
    tassert(Df3[0] == 2*x && Df3[1] == 2*y, "Differentiate dot product");
    tassert(Df4[0] == 1 && Df4[1] == -1,    "Differentiate difference");
    tassert(Df5 == -1,                      "Differentiate negation");
    tassert(approx(Df6,-1/x/x,1e-14),       "Differentiate division");
    tassert(approx(Df7,1,1e-14),            "Differentiate sqrt");
    tassert(approx(Df8,1/x,1e-14),          "Differentiate log");
    tassert(approx(Df9,1,1e-14),            "Differentiate log(exp(x))");
    tassert(approx(Df10,y*exp(x*y),1e-14),  "Differentiate exp");
}


void test_complex()
{
    std::complex<double> z(0,1);
    std::complex<double> zI[4];
    double a, b;
    using std::real;
    using std::imag;

    /* <generator matexpr>
    complex input z;
    output zI(2,2);
    output a, b;
    zI = z*eye(2);
    a = real(z);
    b = imag(z);
    */

    tassert(zI[0] == z && zI[1] == 0.0 && zI[2] == 0.0 && zI[3] == z,
            "Multiply by complex scalar");
    tassert(a == 0 && b == 1,
            "Type on real and imag functions");
}


template<class T>
void test_complexT(T z)
{
    T zI[4];

    /* <generator matexpr complex="T">
    complex input z;
    output zI(2,2);
    zI = z*eye(2);
    */

    tassert(zI[0] == z && zI[1] == 0.0 && zI[2] == 0.0 && zI[3] == z,
            "Multiply by complex scalar");
}


void test_symmetric()
{
    double A[4] = {2, 0,  1, 2};
    double y[2];

    /* <generator matexpr>
    input A symmetric(2);
    output y(2);
    x = [1; 2];
    y = A*x;
    */
    tassert(y[0] == 4 && y[1] == 5, "Load symmetric matrix");
}


void test_funcall()
{
    double a = 1;
    double b;
    /* <generator matexpr>
    input a;
    output b;
    function foo(x) = 10*x;
    b = foo(a) + foo(2);
    */
    tassert(b == 30, "Function call");
}


int main()
{
    test_scalar_mat();
    test_mat_mat();
    test_mat();
    test_opt();
    test_call();
    test_subscript();
    test_inout();
    test_iconst();
    test_range();
    test_deriv();
    test_complex();
    test_complexT(std::complex<double>(0,1));
    test_complexT(double(1));
    test_symmetric();
    test_funcall();
    return test_status;    
}
