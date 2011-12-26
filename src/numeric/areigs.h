#ifndef AREIGS_H
#define AREIGS_H

#include "qcomplex.h"

class Arpack {
public:
    Arpack();
    virtual ~Arpack();

    void set_n(int n_);
    void set_nev(int nev_);
    void set_ncv(int ncv_);
    void set_rvec(int rvec_);
    void set_maxitr(int maxitr_);
    void set_which(char* which_);
    void set_tol(double tol_);

    int    get_n()      { return n;      }
    int    get_nev()    { return nev;    }
    int    get_ncv()    { return ncv;    }
    int    get_rvec()   { return rvec;   }
    int    get_maxitr() { return maxitr; }
    char*  get_which()  { return which;  }
    double get_tol()    { return tol;    }

protected:
    int n;
    int nev;
    int ncv;
    int rvec;
    int mode;
    int maxitr;
    char bmat[2];
    char which[3];
    double tol;
};


class ArpackDS : public Arpack {
public:
    ArpackDS();
    ~ArpackDS();
    int compute_eigs(double* d, double* v);
    void set_mode(int mode_);
    void set_shift(double sigma_);
protected:
    double sigma;
    virtual void times_OP1(double* x, double* opx);
    virtual void times_OP2(double* x, double* opx, double* Mx);
    virtual void times_M  (double* x, double* Mx);
};


class ArpackDN : public Arpack {
public:
    ArpackDN();
    ~ArpackDN();
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
};


class ArpackZN : public Arpack {
public:
    ArpackZN();
    ~ArpackZN();
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
};


#endif /* AREIGS_H */
