/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#ifndef SHAPES_H
#define SHAPES_H

/**
 * @file shapes.h
 *
 * Functions associated with shape function evaluation and
 * isoparametric and PML coordinate transformations.
 */
#include <vector>

#include "qcomplex.h"

class Mesh;


class Shape {
 public:
    Shape(int nen, int ndm);
    Shape();

    virtual ~Shape();
    virtual void eval(double* X) = 0;

    double* N()                { return N_;              }
    double* dNdx()             { return dNdx_;           }
    double& N(int i)           { return N_[i];           }
    double& dNdx(int i, int j) { return dNdx_[i+ndm_*j]; }
    double  x(int i)           { return x_[i];           }
    double  J()                { return J_;              }
    int     nen()              { return nen_;            }
    int     ndm()              { return ndm_;            }

    double   interp(double*   f);
    dcomplex interp(dcomplex* f);
    void     interp(double*   f, double*   fx, int nfields = 1);
    void     interp(dcomplex* f, dcomplex* fx, int nfields = 1);

 protected:
    int owner;
    int ndm_;
    int nen_;
    double  x_[3];
    double* N_;
    double* dNdx_;
    double  J_;
};


class PMLShape : public Shape {
public:
    PMLShape(Shape& base, dcomplex* stretch = NULL);

    dcomplex* dNdxt()             { return &dNdxt_[0];       }
    dcomplex& dNdxt(int i, int j) { return dNdxt_[i+ndm_*j]; }
    dcomplex  Jt()                { return Jt_;              }

    void eval(double* X);

 private:
    Shape& base_;
    std::vector<dcomplex> stretch_;
    std::vector<dcomplex> dNdxt_;
    dcomplex  Jt_;
};


class ShapeInt {
 public:
    virtual ~ShapeInt();
    virtual int end() = 0;
    virtual int next() = 0;

    double* X()      { return X_;    }
    double  X(int i) { return X_[i]; }
    double  W()      { return W_;    }

 protected:
    double X_[3];
    double W_;
};


class Shape1d : public Shape {
 public:
    Shape1d(int nen);
    double& nx(int i) { return nodex[i]; }
    void compute_jacobian(double* dNdX, double* Jmat, double& J);
 protected:
    void remap_gradients(double* Jmat, double* dNdX);
    std::vector<double> nodex;
};


class Shape2d : public Shape {
 public:
    Shape2d(int nen);
    double& nx(int i, int j) { return nodex[2*j+i]; }
    void compute_jacobian(double* dNdX, double* Jmat, double& J);
 protected:
    void remap_gradients(double* Jmat, double* dNdX);
    std::vector<double> nodex;
};


class Shape3d : public Shape {
 public:
    Shape3d(int nen);
    double& nx(int i, int j) { return nodex[3*j+i]; }
    void compute_jacobian(double* dNdX, double* Jmat, double& J);
 protected:
    void remap_gradients(double* Jmat, double* dNdX);
    std::vector<double> nodex;
};


class Quad1d : public Shape1d {
public:
    Quad1d(Mesh* mesh, int eltid, int nen);
    void eval(double* X);
    void inv(double* x, double* X, int maxiter=6, double tol=1e-15);
private:
    void apply_inv(double* J1, double* x, double* xi);
};


class Quad2d : public Shape2d {
public:
    Quad2d(Mesh* mesh, int eltid, int nen);
    void eval(double* X);
    void inv(double* x, double* X, int maxiter=6, double tol=1e-15);
private:
    void apply_inv(double* J1, double* x, double* xi);
};


class Quad3d : public Shape3d {
public:
    Quad3d(Mesh* mesh, int eltid, int nen);
    void eval(double* X);
    void inv(double* x, double* X, int maxiter=6, double tol=1e-15);
private:
    void apply_inv(double* J1, double* x, double* xi);
};


class Quad1dInt : public ShapeInt {
 public:
    Quad1dInt(int nen);
    int end();
    int next();
 private:
    int ix;
    int ngauss;
    void eval();
};


class Quad2dInt : public ShapeInt {
 public:
    Quad2dInt(int nen);
    int end();
    int next();
 private:
    int ix, iy;
    int ngauss;
    void eval();
};


class Quad3dInt : public ShapeInt {
 public:
    Quad3dInt(int nen);
    int end();
    int next();
 private:
    int ix, iy, iz;
    int ngauss;
    void eval();
};


/** Evaluate 1D shape functions and parent gradients
 *
 * @param order    Order of interpolation (nen-1)
 * @param N        Output array of shape functions at x (nen)
 * @param dNdX     Output array of shape gradients at x (nen)
 * @param x        Evaluation point in parent coords
 */
void shape1d(int order, double* N, double* dNdX, double x);


/** Evaluate 2D shape functions and parent gradients
 *
 * @param nen      Number of element nodes
 * @param N        Output array of shape functions at x (nen)
 * @param dNdX     Output array of shape gradients at x (ndm-by-nen)
 * @param x        Evaluation point in parent coords (ndm)
 */
void shape2d(int nen, double* N, double* dNdX, double* x);


/** Evaluate 3D shape functions and parent gradients
 *
 * @param nen      Number of element nodes
 * @param N        Output array of shape functions at x (nen)
 * @param dNdX     Output array of shape gradients at x (ndm-by-nen)
 * @param x        Evaluation point in parent coords (ndm)
 */
void shape3d(int nen, double* N, double* dNdX, double* x);


/** Return order of accuracy for nen-node 2D element
 *
 * @param nen      Number of element nodes
 */
inline int order2d(int nen)
{
    if (nen == 4)  return 1;
    if (nen == 9)  return 2;
    if (nen == 16) return 3;
    if (nen == 25) return 4;
    return 0;
}

/** Return order of accuracy for nen-node 3D element
 *
 * @param nen      Number of element nodes
 */
inline int order3d(int nen)
{
    if (nen == 8)   return 1;
    if (nen == 27)  return 2;
    if (nen == 64)  return 3;
    if (nen == 125) return 4;
    return 0;
}

#endif /* SHAPES_H */
