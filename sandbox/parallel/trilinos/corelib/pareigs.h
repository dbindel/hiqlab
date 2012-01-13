#ifndef PAREIGS_H
#define PAREIGS_H

#include "areigs.h"
#include "qcomplex.h"

class PArpackDS : public Arpack {

 public:

    PArpackDS();
    ~PArpackDS();
    int compute_eigs(double* d, double* v);
    void set_mode(int mode_);
    void set_shift(double sigma_);

 protected:

    double sigma;
    virtual void times_OP1(double* x, double* opx);
    virtual void times_OP2(double* x, double* opx, double* Mx);
    virtual void times_M  (double* x, double* Mx);

    virtual void broadcast(int* x, int n, int root) = 0;
    virtual void broadcast(double* x, int n, int root) = 0;
    virtual void barrier() = 0;
    virtual int  mypid() = 0;
};

class PArpackDN : public Arpack {

 public:

    PArpackDN();
    ~PArpackDN();
    int compute_eigs(double* dr, double* di, double* v);
    void set_mode(int mode_);
    void set_shift(double sigmar_, double sigmai_ = 0);
    void set_shift(dcomplex sigma_) {
        set_shift(real(sigma_), imag(sigma_));
    }

 protected:

    double sigmar, sigmai;
    virtual void times_OP1(double* x, double* opx);
    virtual void times_OP2(double* x, double* opx, double* Mx);
    virtual void times_M  (double* x, double* Mx);

    virtual void broadcast(int* x, int n, int root) = 0;
    virtual void broadcast(double* x, int n, int root) = 0;
    virtual void barrier() = 0;
    virtual int  mypid() = 0;
};

class PArpackZN : public Arpack {

 public:

    PArpackZN();
    ~PArpackZN();
    int compute_eigs(dcomplex* d, dcomplex* v);
    void set_mode(int mode_);
    void set_shift(dcomplex sigma_);
    void set_shift(double sigmar_, double sigmai_ = 0) {
        set_shift(dcomplex(sigmar_, sigmai_));
    }

 protected:

    dcomplex sigma;
    virtual void times_OP1(dcomplex* x, dcomplex* opx);
    virtual void times_OP2(dcomplex* x, dcomplex* opx, dcomplex* Mx);
    virtual void times_M  (dcomplex* x, dcomplex* Mx);

    virtual void broadcast(int* x, int n, int root) = 0;
    virtual void broadcast(double* x, int n, int root) = 0;
    virtual void broadcast(dcomplex* x, int n, int root) = 0;
    virtual void barrier() = 0;
    virtual int  mypid() = 0;
};

#endif /* PAREIGS_H */
