// uses NTL
//   http://www.shoup.net/ntl
// reference: H. Cohen,
//  "A Course in Computational Algebraic Number Theory"
//   section 10.1

#include<NTL/ZZ_p.h>
#include<NTL/mat_GF2.h>
#include "hash_list.h"
using namespace NTL;

#define MB_BOUND 0.5
#define MB_EXTRA 5

struct mb_relation {
    static long len;
    static inline void init(long l) { len=l; } 
    ZZ_p u;
    long *e,t;
    mb_relation() { e = new long[len]; }
    mb_relation(const mb_relation& r) : t(r.t), u(r.u) {
        e = new long[len];
        for(long i=0; i<len; i++) e[i] = r.e[i];
    }
    ~mb_relation() { delete[] e; }
    inline bool operator==(const mb_relation& r) const { return t==r.t; }
    inline bool operator>(const mb_relation& r) const { return t>r.t; }
    inline void operator/=(const mb_relation& r) {
        for(long i=0; i<len; i++) e[i] -= r.e[i];
    }
};

long mb_relation::len;
inline void conv(long& v, const mb_relation& r) { v=r.t; }

void morbri(ZZ& d, const ZZ& n)
// d = factor of n
// by continued fraction method of Morrison and Brillhart
{
    long Jacobi(long, long);
    double log_n = log(conv<double>(n));
    double log_B = MB_BOUND*sqrt(log_n*log(log_n));
    long B = long(exp(log_B));
    long K,Q,i,j,k,p;
    Vec<long> F;
    PrimeSeq ps;
    F.append(ps.next());
    while((p = ps.next())<=B) {
        if((j = Jacobi(n%p,p)) > 0) F.append(p);
        else if(j==0) { d=p; return; }
    }
    K = F.length();
    long N(K+MB_EXTRA), e[K];
    ZZ_p::init(n);
    ZZ_p u(1),v(0),w,x;
    ZZ a(1),b(0),c(n),q,r,s;
    power(q, F[K-1], 2);
    if(q.SinglePrecision()) conv(Q,q);
    else Q = NTL_SP_BOUND;
    mb_relation::init(K+1);
    mb_relation R[N], *rp[N];
    SqrRoot(s,n);
    HashList<mb_relation> T(NumBits(B));
    for(i=0, j=1; i<N; j^=1) {
        add(r,b,s);
        DivRem(q,r,r,a);
        sub(r,s,r);
        b-=r;
        MulAddTo(c,b,q);
        b=r;
        swap(a,c);
        conv(x,q);
        w=u; u*=x; u+=v; v=w;
        for(k=0; k<K; k++) R[i].e[k] = 0;
        q=a;
        for(k=0; k<K; k++) {
            if(F[k]>q) break;
            while(divide(q,q,F[k])) R[i].e[k]++;
        }
        if(q>=Q) continue;
        conv(R[i].t, q);
        R[i].u = u;
        R[i].e[K] = j;
        if(R[i].t == 1) rp[i++] = 0;
        else if(rp[i] = T.lookup_or_install(R[i])) {
            R[i] /= *rp[i];
            i++;
        }
    }
    Mat<GF2> A,X;
    A.SetDims(N,K+1);
    for(i=0; i<N; i++)
        for(j=0; j<=K; j++)
            if(R[i].e[j] & 1) set(A[i][j]);
    kernel(X,A);
    for(k=0; k<X.NumRows(); k++) {
        set(u); set(v);
        for(i=0; i<K; i++) e[i] = 0;
        for(i=0; i<N; i++) {
            if(IsZero(X[k][i])) continue;
            u *= R[i].u;
            if(rp[i]) v *= rp[i]->u;
            for(j=0; j<K; j++) e[j] += R[i].e[j];
        }
        for(i=0; i<K; i++) {
            if(e[i] == 0) continue;
            conv(w, F[i]);
            j = abs(e[i])>>1;
            power(w,w,j);
            if(e[i]>0) v*=w; else u*=w;
        }
        conv(a,u);
        conv(b,v);
        a-=b;
        GCD(d,a,n);
        if(d>1 && d<n) return;
    }
    //	std::cout << "failed to find a factor\n";
    ps.reset(3);
    do {	p = ps.next();
        mul(s,n,p);
        morbri(d,s);
        divide(d,d,p);
    } while(IsOne(d) || d==n);
}
