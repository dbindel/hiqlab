/* HiQLab
 * Copyright (c): Regents of the University of California
 */

#include <cstdio>
#include <cassert>
#include <algorithm>

#include "cscmatrix.h"
#include "qcomplex.h"

#define ME CSCMatrix

using std::fill;
using std::copy;
using std::vector;


ME::ME(int real_flag) :
    m(0), n(0), real_flag(real_flag)
{
}


ME::ME(int m, int n, int nnz, int real_flag) :
    m(m), n(n), real_flag(real_flag),
    jc(n+1), ir(nnz), 
    pr(nnz), pi(real_flag ? 0 : nnz)
{
    fill(jc.begin(), jc.end(), 0);
    fill(ir.begin(), ir.end(), 0);
    fill(pr.begin(), pr.end(), 0);
    if (real_flag)
        fill(pi.begin(), pi.end(), 0);
    jc[n] = nnz;
}


ME::ME(const ME& matrix) :
    m(matrix.m), n(matrix.n), 
    real_flag(matrix.real_flag),
    jc(matrix.jc), ir(matrix.ir),
    pr(matrix.pr), pi(matrix.pi)
{
}


void ME::set_dims(int xm, int xn)
{
    m = xm;
    n = xn;
    jc.resize(n+1);
}


void ME::set_nnz(int xnnz)
{
    jc[n] = xnnz;
    ir.resize(xnnz);
    pr.resize(xnnz);
    if (!real_flag)
        pi.resize(xnnz);
}


void ME::scale(double* du, double* df)
{
    for (int j = 0; j < n; ++j) {
        double dcol = du ? du[j] : 1;
        for (int ii = jc[j]; ii < jc[j+1]; ++ii) {
            double drow = df ? 1.0/df[ir[ii]] : 1;
            pr[ii] *= (dcol * drow);
            if (real_flag)
                pi[ii] *= (dcol * drow);
        }
    }
}


void ME::apply(double* xx, double* xz, double* axx, double* axz)
{
    if (real_flag) {

        // Apply the real operator
        apply(xx, axx);
        apply(xz, axz);

    } else {

        // Apply the complex operator
        fill(axx, axx+m, 0);
        fill(axz, axz+m, 0);
        for (int j = 0; j < n; ++j) {
            double xxj = xx[j];
            double xzj = xz[j];
            for (int i = jc[j]; i < jc[j+1]; ++i) {
                int ii = ir[i];
                axx[ii] += pr[i]*xxj - pi[i]*xzj;
                axz[ii] += pr[i]*xzj + pi[i]*xxj;
            }
        }

    }
}


void ME::apply(dcomplex* x, dcomplex* ax)
{
    fill(ax, ax+m, 0);
    if (real_flag) {

        // Apply the real operator
        for (int j = 0; j < n; ++j) {
            dcomplex xj = x[j];
            for (int i = jc[j]; i < jc[j+1]; ++i) {
                int ii = ir[i];
                ax[ii] += pr[i] * xj;
            }
        }

    } else {

        // Apply the complex operator
        for (int j = 0; j < n; ++j) {
            dcomplex xj = x[j];
            for (int i = jc[j]; i < jc[j+1]; ++i) {
                int ii = ir[i];
                ax[ii] += dcomplex(pr[i], pi[i]) * xj;
            }
        }

    }
}


void ME::apply(double* x, double* ax)
{
    // Real operator only
    assert(real_flag);
    fill(ax, ax+n, 0);
    for (int j = 0; j < n; ++j) {
        double xj = x[j];
        for (int i = jc[j]; i < jc[j+1]; ++i) {
            int ii = ir[i];
            ax[ii] += pr[i] * xj;
        }
    }
}


void ME::dump(const char* fname)
{
    FILE* fp = fopen(fname, "w+");
    assert(fp); // FIXME -- good to fail more gracefully here
    for (int j = 0; j < n; ++j) {
        for (int idx = jc[j]; idx < jc[j+1]; ++idx) {
            fprintf(fp, "%d %d %17.17e", ir[idx]+1, j+1, pr[idx]);
            if (!real_flag)
                fprintf(fp, " %17.17e", pi[idx]);
            fprintf(fp, "\n");
        }
    }
    fclose(fp);
}


void ME::load(const char* fname)
{
    FILE* fp = fopen(fname, "r");
    assert(fp); // FIXME -- good to fail more gracefully here

    fill(jc.begin(), jc.end(), 0);
    ir.resize(0);
    pr.resize(0);
    pi.resize(0);

    int nnz = 0;
    char buf[256];
    while (fgets(buf, sizeof(buf), fp)) {
        int i, j;
        double ar = 0, ai = 0;
        if (sscanf(buf, "%d %d %lf %lf", &i, &j, &ar, &ai) >= 3) {
            jc[j] = ++nnz;
            ir.push_back(i-1);
            pr.push_back(ar);
            if (!real_flag)
                pi.push_back(ai);
        }
    }
    fclose(fp);

    for (int j = 1; j < n+1; ++j)
        if (jc[j] == 0)
            jc[j] = jc[j-1];
}


void ME::load(int* ii, int* jj, double* rdata, double* idata, int nnz)
{
    fill(jc.begin(), jc.end(), 0);
    ir.resize(nnz);
    pr.resize(nnz);
    pi.resize(nnz);

    for (int k = 0; k < nnz; ++k) {
        jc[jj[k]] = k+1;
        ir[k] = ii[k]-1;
        pr[k] = rdata[k];
        if (!real_flag)
            pi[k] = (idata ? idata[k] : 0);
    }

    for (int j = 1; j < n+1; ++j)
        if (jc[j] == 0)
            jc[j] = jc[j-1];
}


/*************************************************************/


/*
 * An assembler for a sparse vector index keeps two data structures: a
 * bit mask that indicates the nonzero pattern, and an index set that
 * indicates where the nonzeros are.  We keep track of the bit mask
 * storage, but the user provides the storage for the index structure.
 * We can do three things with the assembler: start assembling, add an
 * entries, and clear.
 */

class SpVecStructAssembler {
public:
    SpVecStructAssembler(int m) : mark_(m), ir_(NULL), nnz_(0) {}

    void start(vector<int>& ir, int start = 0) {
        clear();
        ir_ = &(ir[start]);
    }

    void add(vector<int>& irx, int l, int u) {
        for (int ii = l; ii < u; ++ii) {
            int i = irx[ii];
            if (!mark_[i]) {
                mark_[i] = true;
                ir_[nnz_++] = i;
            }
        }
    }

    void sort() {
        std::sort(ir_, ir_+nnz_);
    }

    void clear() { 
        if (ir_ != NULL)
            for (int ii = 0; ii < nnz_; ++ii)
                mark_[ir_[ii]] = 0;
        nnz_ = 0;
    }

    int nnz() { return nnz_; }

private:
    vector<bool> mark_;  // Bit vector representing sparsity
    int*         ir_;    // Indices of marked nonzeros
    int          nnz_;   // Number of nonzeros
};


/*
 * An assembler for a sparse vector keeps a full representation of the
 * vector, but also has a (user-provided) list of nonzero indices to
 * make copying and clearing the entries faster.
 */

class SpVecAssembler {
public:
    SpVecAssembler(int m) : data_(m), ir_(NULL) {}

    void start(vector<int>& ir, int l, int u) {
        if (ir_)
            clear();
        ir_ = &(ir[0]);
        l_ = l;
        u_ = u;
    }

    void finish(vector<double>& pr) {
        for (int i = l_; i < u_; ++i)
            pr[i] = data_[ir_[i]];
    }

    void add(vector<int>& ir, vector<double>& pr, double alpha, int l, int u) {
        for (int i = l; i < u; ++i)
            data_[ir[i]] += alpha*pr[i];
    }

    void clear() {
        for (int i = l_; i < u_; ++i)
            data_[ir_[i]] = 0;
    }

    double operator[](int i) { return data_[i]; }

private:
    vector<double> data_;    // Full vector storage
    int*           ir_;      // Nonzeros are in ir_[l_:u_-1]
    int            l_, u_;
};


int ME::mul_col_struct(SpVecStructAssembler& Axs, ME& X, int k)
{
    for (int jj = X.jc[k]; jj < X.jc[k+1]; ++jj) {
        int j = X.ir[jj];
        Axs.add(ir, jc[j], jc[j+1]);
    }
    return Axs.nnz();
}


int ME::mul_col_struct(SpVecStructAssembler& Axs, vector<int>& irx, int nnzx)
{
    for (int jj = 0; jj < nnzx; ++jj) {
        int j = irx[jj];
        Axs.add(ir, jc[j], jc[j+1]);
    }
    return Axs.nnz();
}


void ME::mul_col_real(SpVecAssembler& Axr, ME& X, int k)
{
    for (int jj = X.jc[k]; jj < X.jc[k+1]; ++jj) {
        int j = X.ir[jj];
        Axr.add(ir, pr, X.pr[jj], jc[j], jc[j+1]);
        if (!real_flag && !X.real_flag)
            Axr.add(ir, pi, -X.pi[jj], jc[j], jc[j+1]);
    }
}


void ME::mul_col_real(SpVecAssembler& Axr, 
                      SpVecAssembler& xr, SpVecAssembler& xi,
                      bool real_x, vector<int>& irx, int nnzx)
{
    for (int jj = 0; jj < nnzx; ++jj) {
        int j = irx[jj];
        Axr.add(ir, pr, xr[j], jc[j], jc[j+1]);
        if (!real_flag && !real_x)
            Axr.add(ir, pi, -xi[j], jc[j], jc[j+1]);
    }
}


void ME::mul_col_imag(SpVecAssembler& Axi, ME& X, int k)
{
    for (int jj = X.jc[k]; jj < X.jc[k+1]; ++jj) {
        int j = X.ir[jj];
        if (!real_flag)   Axi.add(ir, pi, X.pr[jj], jc[j], jc[j+1]);
        if (!X.real_flag) Axi.add(ir, pr, X.pi[jj], jc[j], jc[j+1]);
    }
}


void ME::mul_col_imag(SpVecAssembler& Axr, 
                      SpVecAssembler& xr, SpVecAssembler& xi,
                      bool real_x, vector<int>& irx, int nnzx)
{
    for (int jj = 0; jj < nnzx; ++jj) {
        int j = irx[jj];
        if (!real_flag) Axr.add(ir, pi, xr[j], jc[j], jc[j+1]);
        if (!real_x)    Axr.add(ir, pr, xi[j], jc[j], jc[j+1]);
    }
}


int ME::nnz_mul(ME& X)
{
    SpVecStructAssembler Axs(m);
    vector<int> irAxs(m);
    int nnz = 0;
    for (int k = 0; k < X.n; ++k) {
        Axs.start(irAxs);
        nnz += mul_col_struct(Axs, X, k);
    }
    return nnz;
}


int ME::nnz_mul(ME& Y, ME& X)
{
    SpVecStructAssembler Axs(m);
    SpVecStructAssembler YAxs(Y.m);
    vector<int> irAx(m);
    vector<int> irYAx(Y.m);
    int nnz = 0;
    for (int k = 0; k < X.n; ++k) {
        Axs.start(irAx);
        YAxs.start(irYAx);
        nnz += mul_col_struct(YAxs, irAx, mul_col_struct(Axs, X, k));
    }
    return nnz;
}


void ME::mul(ME& X, ME& AX, bool build_index)
{
    SpVecStructAssembler Axs(m);
    SpVecAssembler Ax(m);
    AX.jc[0] = 0;
    for (int k = 0; k < X.n; ++k) {
        if (build_index) {
            Axs.start(AX.ir, AX.jc[k]);
            AX.jc[k+1] = AX.jc[k] + mul_col_struct(Axs, X, k);
            Axs.sort();
        } 
            
        Ax.start(AX.ir, AX.jc[k], AX.jc[k+1]);
        mul_col_real(Ax, X, k);
        Ax.finish(AX.pr);

        if (!AX.real_flag) {
            Ax.start(AX.ir, AX.jc[k], AX.jc[k+1]);
            mul_col_imag(Ax, X, k);
            Ax.finish(AX.pi);
        }
    }
}


void ME::mul(ME& Y, ME& X, ME& YAX, bool build_index)
{
    SpVecStructAssembler Axs(m);
    SpVecAssembler Axr(m);
    SpVecAssembler Axi(m);
    vector<int> irAx(m);
    bool real_Ax = (real_flag && X.real_flag);

    SpVecStructAssembler YAxs(Y.m);
    SpVecAssembler YAx(m);

    YAX.jc[0] = 0;
    for (int k = 0; k < X.n; ++k) {
        Axs.start(irAx);
        int nnzAx = mul_col_struct(Axs, X, k);
        Axr.start(irAx, 0, nnzAx);
        mul_col_real(Axr, X, k);
        if (real_Ax) {
            Axi.start(irAx, 0, nnzAx);
            mul_col_imag(Axi, X, k);
        }

        if (build_index) {
            YAxs.start(YAX.ir, YAX.jc[k]);
            YAX.jc[k+1] = YAX.jc[k] + Y.mul_col_struct(YAxs, irAx, nnzAx);
            YAxs.sort();
        }

        YAx.start(YAX.ir, YAX.jc[k], YAX.jc[k+1]);
        Y.mul_col_real(YAx, Axr, Axi, real_Ax, irAx, nnzAx);
        YAx.finish(YAX.pr);

        if (!YAX.real_flag) {
            YAx.start(YAX.ir, YAX.jc[k], YAX.jc[k+1]);
            Y.mul_col_imag(YAx, Axr, Axi, real_Ax, irAx, nnzAx);
            YAx.finish(YAX.pi);
        }

    }
}


// Compute matmul, allocating a new matrix if needed
ME* ME::mul(ME* B, ME* result)
{
    if (!result) {
        int nnz = nnz_mul(*B);
        result = new ME(m, B->n, nnz, real_flag && B->real_flag);
        mul(*B, *result, true);
    } else {
        assert(n == B->m);
        assert(!result->real_flag || (real_flag && B->real_flag));
        mul(*B, *result, false);
    }
    return result;
}


// Compute triple product
ME* ME::mul(ME* Y, ME* X, ME* result)
{
    assert(Y->n == m && n == X->m);
    if (!result) {
        int nnz = nnz_mul(*Y, *X);
        result = new ME(Y->m, X->n, nnz, 
                        real_flag && X->real_flag && Y->real_flag);
        mul(*Y, *X, *result, true);
    } else {
        assert(result->m == m && result->n == X->n);
        assert(!result->real_flag || 
               (real_flag && X->real_flag && Y->real_flag));
        mul(*Y, *X, *result, false);
    }
    return result;
}


int ME::nnz_add(ME& X)
{
    SpVecStructAssembler sums(m);
    vector<int> irs(m);
    int nnz = 0;
    for (int k = 0; k < X.n; ++k) {
        sums.start(irs);
        sums.add(ir, jc[k], jc[k+1]);
        sums.add(X.ir, X.jc[k], X.jc[k+1]);
        nnz += sums.nnz();
    }
    return nnz;
}


void ME::add(double alpha, double beta, ME& B, ME& C, bool build_index)
{
    SpVecStructAssembler cs(m);
    SpVecAssembler cx(m);
    C.jc[0] = 0;
    for (int k = 0; k < n; ++k) {
        if (build_index) {
            cs.start(C.ir, C.jc[k]);
            cs.add(ir, jc[k], jc[k+1]);
            cs.add(B.ir, B.jc[k], B.jc[k+1]);
            C.jc[k+1] = C.jc[k] + cs.nnz();
            cs.sort();
        } 

        cx.start(C.ir, C.jc[k], C.jc[k+1]);
        cx.add(ir,   pr,   alpha, jc[k],   jc[k+1]);
        cx.add(B.ir, B.pr, beta,  B.jc[k], B.jc[k+1]);
        cx.finish(C.pr);

        if (!C.real_flag) {
            cx.start(C.ir, C.jc[k], C.jc[k+1]);
            if (!real_flag)   cx.add(ir,   pi,   alpha, jc[k],   jc[k+1]);
            if (!B.real_flag) cx.add(B.ir, B.pi, beta,  B.jc[k], B.jc[k+1]);
            cx.finish(C.pi);
        }
    }
}


ME* ME::add(double alpha, double beta, ME* B, ME* result)
{
    assert(n == B->n && m == B->m);
    if (!result) {
        int nnz = nnz_add(*B);
        result = new ME(m, n, nnz, real_flag && B->real_flag);
        add(alpha, beta, *B, *result, true);
    } else {
        assert(n == result->n && m == result->m);
        add(alpha, beta, *B, *result, false);
    }
    return result;
}


/*************************************************************/


ME* ME::transpose()
{
    int nnz = jc[n];
    ME* result = new ME(n, m, nnz, real_flag);

    int* new_jc = result->get_jc();
    int* new_ir = result->get_ir();
    double* new_pr = result->get_Ax();
    double* new_pi = result->get_Az();
    fill(result->jc.begin(), result->jc.end(), 0);

    // Count nonzeroes per column (offset by one)
    for (int k = 0; k < nnz; ++k)
        ++new_jc[ir[k]+1];

    // Accumulate counts
    for (int k = 1; k < m+1; ++k)
        new_jc[k] += new_jc[k-1];

    // Bucket sort data into place
    for (int j = 0; j < n; ++j) {
        for (int ii = jc[j]; ii < jc[j+1]; ++ii) {
            int slot = new_jc[ir[ii]]++;
            new_ir[slot] = j;
            new_pr[slot] = pr[ii];
            if (!real_flag)
                new_pi[slot] = pi[ii];
        }
    }

    // Shift counts
    for (int k = m-1; k >= 0; --k)
        new_jc[k+1] = new_jc[k];
    new_jc[0] = 0;

    return result;
}


ME* ME::submatrix(int i1, int i2, int j1, int j2)
{
    assert(0 <= i1 && i1 <= i2 && i2 <= n);
    assert(0 <= j1 && j1 <= j2 && j2 <= n);

    int new_m = (i2-i1)+1;
    int new_n = (j2-j1)+1;
    ME* result = new ME(new_m, new_n, 0, real_flag);

    // Set up column index structure
    int new_nnz = 0;
    int* new_jc = result->get_jc();
    for (int j = 0; j < new_n; ++j) {
        new_jc[j] = new_nnz;
        for (int i = jc[j1+j]; i < jc[j1+j+1]; ++i)
            if (i1 <= ir[i] && ir[i] <= i2)
                ++new_nnz;
    }
    new_jc[new_n] = new_nnz;
    result->set_nnz(new_nnz);

    // Set up row indices and data arrays
    int* new_ir = result->get_ir();
    double* new_pr = result->get_Ax();
    double* new_pi = result->get_Az();
    int* irp = new_ir;
    double* prp = new_pr;
    double* pip = new_pi;
    for (int j = 0; j < new_n; ++j) {
        for (int i = jc[j1+j]; i < jc[j1+j+1]; ++i) {
            if (i1 <= ir[i] && ir[i] <= i2) {
                *irp++ = ir[i];
                *prp++ = pr[i];
                if (real_flag)
                    *pip++ = pi[i];
            }
        }
    }

    return result;
}


ME* vcat(ME* A, ME* B)
{
    assert(A->n == B->n);
    int Annz = A->jc[A->n];
    int Bnnz = B->jc[B->n];

    int n = A->n;
    int m = A->m + B->m;
    int nnz = Annz + Bnnz;
    ME* result = new ME(m, n, nnz, A->real_flag && B->real_flag);

    int*    jc = result->get_jc();
    int*    ir = result->get_ir();
    double* pr = result->get_Ax();
    double* pi = result->get_Az();

    int*    irp = ir;
    double* prp = pr;
    double* pip = pi;
    for (int j = 0; j < n; ++j) {
        jc[j] = A->jc[j] + B->jc[j];
        for (int i = A->jc[j]; i < A->jc[j+1]; ++i) {
            *irp++ = A->ir[i];
            *prp++ = A->pr[i];
            if (pi)
                *pip++ = (A->real_flag ? A->pi[i] : 0);
        }
        for (int i = B->jc[j]; i < A->jc[j+1]; ++i) {
            *irp++ = B->ir[i] + A->m;
            *prp++ = B->pr[i];
            if (pi)
                *pip++ = (B->real_flag ? B->pi[i] : 0);
        }
    }

    return result;
}


ME* hcat(ME* A, ME* B)
{
    assert(A->m == B->m);
    int Annz = A->jc[A->n];
    int Bnnz = B->jc[B->n];

    int n = A->n + B->n;
    int m = A->m;
    int nnz = Annz + Bnnz;
    ME* result = new ME(m, n, nnz, A->real_flag && B->real_flag);

    int*    jc = result->get_jc();
    int*    ir = result->get_ir();
    double* pr = result->get_Ax();
    double* pi = result->get_Az();

    copy(&(ir[0]),    &(ir[Annz]),      &(A->ir[0]));
    copy(&(pr[0]),    &(pr[Annz]),      &(A->pr[0]));
    copy(&(ir[Annz]), &(ir[Annz+Bnnz]), &(B->ir[0]));
    copy(&(pr[Annz]), &(pr[Annz+Bnnz]), &(B->pr[0]));
    if (pi) {
        if (A->real_flag) copy(&(pi[0]),    &(pi[Annz]),      &(A->pi[0])); 
        else              fill(&(pi[0]),    &(pi[Annz]),      0);
        if (B->real_flag) copy(&(pi[Annz]), &(pi[Annz+Bnnz]), &(B->pi[0]));
        else              fill(&(pi[Annz]), &(pi[Annz+Bnnz]), 0);
    }

    copy(jc,      jc+A->n,   &(A->jc[0]));
    copy(jc+A->n, jc+B->n+1, &(B->jc[0]));
    for (int j = A->n; j <= n; ++j)
        jc[j] += Annz;

    return result;
}
