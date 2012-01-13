#include <vector>

#include "petsc_jdqz_util.h"
#include "qlapack_util.h"

using std::complex;
using std::real;
using std::imag;

/*@t -----
 * \section{Conduct Projected Gram-Schmidt}
 *
 * Given vectors V, a projection removing components of vectors Z,
 * followed by a QR factorization through Gramm-Schmidt is conducted.
 *
 * 1. gamma=0
 *
 *    (I-Z*Z^H)V     = Q 
 *     V - Z*(Z^H*V) = Q
 *     V - Z*R_Z     = Q
 *
 * 2. gamma=1
 *
 *     V - Z*R_Z     = Q*R_Q
 *
 *     This implies also, V = [Z,Q][R_Z;R_Q]; and R_Q is upper triangular.
 *
 * In both cases,
 *     1. if Q is not given, V is overwritten with Q.
 *     2. if Z is not given, no projection is conducted.
 *     3. if R is not given, R is not stored.
 *
 * 1. type=0 -> CGS
 * 2. type=0 -> MGS
 *@c*/


void ProjectedGS(Vec* V, int size_v, Vec* Z, int size_z,
                Vec* Q, PetscScalar* R, int ldR, Vec tempv, 
                int type, int gamma, int maxiter)
{
    // -- Allocate necessary workspace
    Vec* Qt;
    PetscScalar* vals;
    PetscScalar* valu;
    PetscScalar  valu1;
    vals = new PetscScalar[size_v+size_z];
    valu = new PetscScalar[size_v+size_z];

    // -- Qt points to Q or V for computation
    Qt = new Vec[size_v];
    for (int i = 0; i < size_v; ++i)
        if(Q) {
            VecCopy(V[i],Q[i]);
            Qt[i] = Q[i];
        } else
            Qt[i] = V[i];

    for (int i = 0; i < size_v; ++i) {

        Vec q = Qt[i];
        PetscReal nr_o, nr, nr_n;
        VecNorm(q,NORM_2,&nr_o);
        nr  = 1e-16*nr_o;
        nr_n=   0.1*nr_o;

        memset(vals, 0, (size_v+size_z)*sizeof(PetscScalar));

        // -- Conduct iteration while 
        //    1. Norm is still decreasing (nr_n < 0.5*nr_o)
        //    2. Norm is still not relatively close to zero
        //    3. Iteration number is smaller than  maxiter
        int niter = 0;
        while (nr_n < 0.5*nr_o && nr_n > nr && niter < maxiter) {

            // -- Conduct orthogonalization
            if (type==0) { // -- CGS

                // -- Start inner products for Z
                if (Z) {
                    for (int j = 0; j < size_z; ++j)
                        VecDotBegin(q,Z[j],&valu[j]);                    
                    for (int j = 0; j < size_z; ++j)
                        VecDotEnd  (q,Z[j],&valu[j]);                    

                    // -- Update Z contribution
                    VecSet(tempv, 0.0);
                    VecMAXPY(tempv,size_z,valu,Z);
                    VecAXPY(q,-1.0,tempv);
                    for (int j = 0; j < size_z; ++j)
                        vals[j]+=valu[j];
                }

                // -- Start inner products for Q
                if (gamma==1 && i > 0) {
                    for (int j = 0; j < i; ++j)
                        VecDotBegin(q,Qt[j],&valu[size_z+j]);
                    for (int j = 0; j < i; ++j)
                        VecDotEnd  (q,Qt[j],&valu[size_z+j]);                    

                    // -- Update Q contribution
                    VecSet(tempv, 0.0);
                    VecMAXPY(tempv,i,&valu[size_z],Qt);
                    VecAXPY(q,-1.0,tempv);
                    for (int j = 0; j < i; ++j)
                        vals[size_z+j]+=valu[size_z+j];
                }

            } else if (type==1) { // -- MGS

                // -- Start inner products for Z
                if (Z)
                    for (int j = 0; j < size_z; ++j) {
                        VecDot(q,Z[j],&valu1);
                        VecAXPY(q,-valu1,Z[j]);
                        vals[j]+=valu1;
                    }                  

                // -- Start inner products for Q
                if (gamma==1 && i > 0)
                    for (int j = 0; j < i; ++j) {
                        VecDot(q,Qt[j],&valu1);
                        VecAXPY(q,-valu1,Qt[j]);
                        vals[size_z+j]+=valu1;
                    }

            }

            // -- Update norms
            nr_o = nr_n;
            VecNorm(q,NORM_2,&nr_n);
            niter++;
        }


        // -- Update R
        if (gamma==0) {

            if (R)
                for (int j = 0; j < size_z; ++j)
                    R[ldR*i+j] = vals[j];

        } else if (gamma=1) {

            VecScale(q,PetscScalar(1/nr_n)); 
            vals[size_z+i] += PetscScalar(nr_n);

            if (R)
                for (int j = 0; j < size_z+i+1; ++j)
                    R[ldR*i+j] = vals[j];
        }

    }

    // -- Clean up
    delete[] Qt;
    delete[] valu;
    delete[] vals;
}


/*@t ----------
 * \section{Create initial Arnoldi space}
 *
 *@c*/


void InitArnoldi(Vec* V, Vec v0, Mat A, Mat B, PetscScalar* sigma, int minq)
{
    // -- allocate space for projection
    Vec tempv;
    VecDuplicate(v0,&tempv);

    // -- Set the shift
    PetscScalar tsigma[2];
    if (PetscRealPart(sigma[0])==0 && PetscImaginaryPart(sigma[1])==0) {
        tsigma[0] = PetscScalar(0);
        tsigma[1] = PetscScalar(1);
    } else {
        tsigma[0] = sigma[0];
        tsigma[1] = sigma[1];
    }

    // -- Normalize initial vector
    PetscReal normv;
    VecNorm(v0,NORM_2,&normv);
    VecCopy(v0,V[0]);
    VecScale(V[0],PetscScalar(1/normv));

    // -- Conduct GS
    for (int i = 0; i < minq-1; ++i) {

        HarmonicVector(A,B,sigma,V[i],V[i+1],tempv);

        // -- Orthonormalize
        ProjectedGS(&V[i+1], 1, V, i+1, 
                    NULL, NULL, 0, tempv, 1, 1, 3);
    }

    // -- clean up
    VecDestroy(tempv);
}


/*@t ----------
 * \section{Compute the projection of vectors}
 *@c*/


void ProjectedOP(Mat A, Vec* V, int size_v, Vec* W, int size_w, Vec tempv,
                PetscScalar* MA, int ldMA, int imin, int jmin)
{
    // -- MA[ldMA*(jmin+j)+(imin+i)] = W[i]'*A*V[j]
    for (int j = 0; j < size_v; ++j) {

        // -- A*V[j]
        if (A==PETSC_NULL)
            VecCopy(V[j],tempv);
        else
            MatMult(A, V[j], tempv);

        // -- W[i]'*A*V[j] 
        //FIXME: Valgrind preferes unsplit version??
        for (int i = 0; i < size_w; ++i)
            VecDot(tempv,W[i],&MA[ldMA*(jmin+j)+(imin+i)]);
//        for (int i = 0; i < size_w; ++i)
//            VecDotBegin(tempv,W[i],&MA[ldMA*(jmin+j)+(imin+i)]);
//        for (int i = 0; i < size_w; ++i)
//            VecDotEnd  (tempv,W[i],&MA[ldMA*(jmin+j)+(imin+i)]);
    }

}


/*@t ----------
 * \section{Compute harmonic vector}
 * 
 * Given matrices A,B and pair(sigma[0], sigma[1]) where
 * the eigenvalue is sigma[0]/sigma[1] compute the harmonic vector
 *
 * w = sigma[1]*A*v - sigma[0]*B*v
 *
 *@c*/


void HarmonicVector(Mat A, Mat B, PetscScalar* sigma, Vec v, Vec w, Vec Av)
{
    // -- Av = A*v
    MatMult(A,v,Av);

    // -- w = B*v
    MatMult(B,v,w);

    // -- w = sigma[1]*Av  - sigma[0]*w
    //      = sigma[1]*A*v - sigma[0]*B*v 
    VecAXPBY(w,sigma[1],-sigma[0],Av);
}


/*@t ----------
 * \section{Utility functions implementing LAPACK routines}
 *
 *@c*/

/*@t ----------
 * \section{Utility structures and sorters}
 *
 *@c*/


struct val_ind {
    double val;
    int ind;
};


inline static int compare_val_ind(const val_ind& pair1, const val_ind& pair2)
{
    if (pair1.val < pair2.val) return 1;
    return 0;
}


/*@t ----------
 * \section{Scale eigenvalues}
 *
 * Scale so that eigenvalue pairs (alpha,beta) have modulus between 
 * 0 an 1. Additionally the angle is modified so that imag(beta)=0.
 * The eigenvalues are essentially,
 *
 *    eigenvalue = alpha/beta
 *
 *@c*/


void ScaleEig(complex<double>* a, complex<double>* b, int n)
{
    for (int i = 0; i < n; ++i) {

        // -- magab= sqrt(|a|^2 + |b|^2)
        // -- magb = |b|
        double magab = sqrt( real(a[i])*real(a[i]) + imag(a[i])*imag(a[i])
                            +real(b[i])*real(b[i]) + imag(b[i])*imag(b[i]));
        double magb  = sqrt( real(b[i])*real(b[i]) + imag(b[i])*imag(b[i]));
        if (magb==0) magb = 1;

        // -- a*conj(b)
        complex<double> acb = a[i]*conj(b[i]);

        // -- a*conj(b)/|b|/magab
        // -- |b|/magab
        a[i] = acb/complex<double>(magb*magab, 0);
        b[i] = complex<double>(magb/magab,     0);
    }
}


void ScaleEig(double* ar, double* ai, double* br, double* bi, int n)
{
    complex<double>* a = new complex<double>[n];
    complex<double>* b = new complex<double>[n];

    for (int i = 0; i < n; ++i) {
        a[i] = complex<double>(ar[i],ai[i]);
        b[i] = complex<double>(br[i],bi[i]);
    }
    ScaleEig(a,b,n);
    for (int i = 0; i < n; ++i) {
        ar[i] = real(a[i]);
        ai[i] = imag(a[i]);
        br[i] = real(b[i]);
        bi[i] = imag(b[i]);
    }

    // -- cleanup
    delete[] a;
    delete[] b;
}


/*@t ----------
 * \section{Sort eigenvalues}
 *
 * Sort given eigenvalues in the complex plane based on,
 *
 * 1. the absolute norm of the distance from a given point
 * 2. the chordal metric between two complex pairs (s,t) and (a,b) 
 *    is defined by the acute angle tau between the two vectors (s,t) and (a,b)
 *    as
 *         sin(tau) = sqrt( 1 - cos(tau)^2)
 *                  = sqrt( 1 - (normalize([a,b])*normalize([s,t])')^2)
 *
 *    Here we only would only like to sort. The order is still preserved by modifying
 *    the norm to,
 *
 *         sqrt(1- sin(tau)^2) = normalize([a,b])*normalize([s,t])'
 *
 *    or
 *        sqrt(1 - sin(tau)^2) * norm([a,b]) = [a,b]'*normalize([s,t])'
 *
 * Output is given in [index] as the new ordering.
 *@c*/


void SortEig(complex<double>* s, complex<double>* t, 
             complex<double>  a, complex<double>  b,
             int* index, int n, int sortCD)
{
    std::vector<val_ind> pairs;
    val_ind vi;

    // -- [a, b]
    // -- s = [ s[1]; s[2]; ... s[n]; ];
    // -- t = [ t[1]; t[2]; ... t[n]; ];
    if (sortCD) {

        for (int i = 0; i < n; ++i) {

            // -- mag = sqrt(|s[i]|^2 + |t[i]|^2)
            double mag = sqrt( real(s[i])*real(s[i]) + imag(s[i])*imag(s[i])
                             + real(t[i])*real(t[i]) + imag(t[i])*imag(t[i]) );

            // -- satb = s[i]*conj(a) + t[i]*conj(b);
            complex<double> satb = s[i]*conj(a) + t[i]*conj(b);
            double adot = abs(satb);

            vi.val = -adot/mag;
            vi.ind = i;
            pairs.push_back(vi);
        }

    } else {

        // -- adb = a/b
        complex<double> adb = a/b;

        for (int i = 0; i < n; ++i) {

            // -- sdt = s[i]/t[i]
            complex<double> sdt = s[i]/t[i];

            if (real(b)==0 && imag(b)==0) { 

                // -- vi.val = -abs(s[i]/t[i])
                vi.val = -abs(sdt);

            } else {

                // -- vi.val = abs(s[i]/t[i]-a/b)
                vi.val = abs(sdt-adb);
            }
            vi.ind = i;
            pairs.push_back(vi);
        }

    }

    std::sort(pairs.begin(),pairs.end(),compare_val_ind);

    for (int i = 0; i < n; ++i)
        index[pairs[i].ind] = i;


}


void SortEig(double* sr, double* si, double* tr, double* ti,
             double ar, double ai, double br, double bi,
             int* index, int n, int sortCD)
{
    complex<double>* s = new complex<double>[n];
    complex<double>* t = new complex<double>[n];
    complex<double> a(ar,ai);
    complex<double> b(br,bi);

    for (int i = 0; i < n; ++i) {
        s[i] = complex<double>(sr[i],si[i]);
        t[i] = complex<double>(tr[i],ti[i]);
    }

    SortEig(s,t,a,b,index,n,sortCD);

    // -- cleanup
    delete[] s;
    delete[] t;
}


/*@t ----------
 * \section{Sorting the QZ}
 *
 * Sort the QZ obtained so that the first [gamma] eigenvalues closest
 * to the pair (alpha,beta) are moved to the top-left corner. The movement
 * of the eigenvalues depend on SwapQZ. Currently this method moves the eigenvalues
 * but does not order them.
 *@c*/


void SortQZ(int n, complex<double>* Q, int ldQ,
                   complex<double>* Z, int ldZ,
                   complex<double>* S, int ldS,
                   complex<double>* T, int ldT,
                   complex<double>* MA,int ldMA,
                   complex<double>* MB,int ldMB,
                   complex<double> alpha, complex<double> beta,
                   int gamma)
{
    int index[n];

    // -- nsort = max(1, min( n-1, gamma) )
    int nsort = n - 1;
    nsort = (nsort > gamma)  ? gamma : nsort;
    nsort = (nsort > 1     ) ? nsort  : 1;

    // -- [S,T,Z,Q] = qz(MA,MB); Z = Z';
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i) {
            S[j*ldS+i] = MA[j*ldMA+i];
            T[j*ldT+i] = MB[j*ldMB+i];
        }
    complex<double>* a = new complex<double>[n];
    complex<double>* b = new complex<double>[n];
    zgges(n, S, ldS, T, ldT, a, b, Z, ldZ, Q, ldQ);
    delete[] a;
    delete[] b;

    // -- Scale eigenvalues so that diag(T)>=0
    complex<double>* Di   = new complex<double>[n*n];
    complex<double>* tmat = new complex<double>[n*n];
    memset(Di  , 0, n*n*sizeof(complex<double>));
    memset(tmat, 0, n*n*sizeof(complex<double>));
    for (int i = 0; i < n; ++i) {
        double rT   = real(T[i*ldT+i]);
        double iT   = imag(T[i*ldT+i]);
        double magT = sqrt(rT*rT + iT*iT);
        Di[n*i+i] = complex<double>(rT/magT,-iT/magT);
    }
    zgemm("NN", n, n, n, complex<double>(1,0), Q, ldQ, Di, n,
                         complex<double>(0,0), tmat, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            Q[j*ldQ+i] = tmat[j*n+i];

    zgemm("NN", n, n, n, complex<double>(1,0), S, ldS, Di, n,
                         complex<double>(0,0), tmat, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            S[j*ldS+i] = tmat[j*n+i];

    zgemm("NN", n, n, n, complex<double>(1,0), T, ldT, Di, n,
                         complex<double>(0,0), tmat, n);
    for (int j = 0; j < n; ++j)
        for (int i = 0; i < n; ++i)
            T[j*ldT+i] = tmat[j*n+i];

    delete[] Di;
    delete[] tmat;


    // -- Sort eigenvalues
    complex<double>* s = new complex<double>[n];
    complex<double>* t = new complex<double>[n];
    for (int i = 0; i < n; ++i) {
        s[i] = S[i*ldS+i];
        t[i] = T[i*ldT+i];
    }
    int sortCD = 1;
    SortEig(s,t,alpha,beta,index,n,sortCD);
    delete[] s;
    delete[] t;

    // -- Swap the QZ
    SwapQZ(n,Q,ldQ,Z,ldZ,S,ldS,T,ldT,index,nsort);
}


void SwapQZ(int n, complex<double>* Q, int ldQ,
                   complex<double>* Z, int ldZ,
                   complex<double>* S, int ldS,
                   complex<double>* T, int ldT,
                   int* index, int nsort)
{
    int* select = new int[n];
    memset(select, 0, n*sizeof(int));
    
    // -- Move first 'nsort' to the front 
    for (int i = 0; i < n; ++i)
        if (index[i]<nsort)
            select[i] = 1;

    complex<double>* a = new complex<double>[n];
    complex<double>* b = new complex<double>[n];
    ztgsen(n, S, ldS, T, ldT, Z, ldZ, Q, ldQ, a, b, select);
    delete[] a;
    delete[] b;

    // -- cleanup
    delete[] select;
}
