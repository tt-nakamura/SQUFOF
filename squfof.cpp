// uses NTL
//   http://www.shoup.net/ntl
// reference: H. Cohen,
//  "A Course in Computational Algebraic Number Theory
//   section 8.7

#include<NTL/ZZ.h>
#include "hash_tab.h"
using namespace NTL;

#define SQUFOF_MAX_BITS 63

long squfof(const ZZ& n)
// return a factor of n
// by Square Form Factorization method of Daniel Shanks
// assume n is composite
{
    long IsSquare(ZZ&, const ZZ&);
    long IsSquare(long);
    if(!IsOdd(n)) return 2;
    ZZ D(n),d;
    if(bit(n,1)) D<<=2;
    SqrRoot(d,D);
    if(NumBits(d) > SQUFOF_MAX_BITS)
        TerminalError("squfof: n is too large");
    ldiv_t qr;
    long a,b,c,s,t,L,M;
    long &q(qr.quot), &r(qr.rem);
    conv(s,d);
    if(IsSquare(d,n)) return conv<long>(d);
    a=1; b=s;
    if(IsOdd(D)^(b&1)) b--;
    conv(d,b); d*=b; sub(d,D,d); d>>=2;
    conv(c,d);
    // D=b^2+4ac, a>0,b>0,c>0
    L = SqrRoot(s);
    HashTable<long> Q(NumBits(L));
    if(b&1) {
        s++;
        s>>=1; t=s-1;
        b>>=1; r=b;
        for(;;) {
            r+=s;
            qr = ldiv(r,c);
            r=t-r;
            b-=r;
            b*=q;
            a+=b;
            b=r;
            if(a==1) goto a;
            if(c<=L) Q.install(c);
            r+=s;
            qr = ldiv(r,a);
            r=t-r;
            b-=r;
            b*=q;
            c+=b;
            b=r;
            if((q=IsSquare(a)) && Q.lookup(q)==0) break;
            if(c==1) goto a;
            if(a<=L) Q.install(a);
        }
        if((a = GCD(q, (b<<1)|1)) > 1) return a;
        for(a=q, c*=q;;) {
            r+=s;
            qr = ldiv(r,a);
            r=t-r; if(r==b) return a;
            b-=r;
            b*=q;
            c+=b;
            b=r;
            r+=s;
            qr = ldiv(r,c);
            r=t-r; if(r==b) return c;
            b-=r;
            b*=q;
            a+=b;
            b=r;
        }
    }
    else {
        M=(L<<1);
        s>>=1;
        b>>=1; r=b;
        for(;;) {
            r+=s;
            qr = ldiv(r,c);
            r=s-r;
            b-=r;
            b*=q;
            a+=b;
            b=r;
            if(a==1) goto a;
            if(c<=M) {
                if((c&1)==0) Q.install(c>>1);
                else if(c<=L) Q.install(c);
            }
            r+=s;
            qr = ldiv(r,a);
            r=s-r;
            b-=r;
            b*=q;
            c+=b;
            b=r;
            if((q=IsSquare(a)) && Q.lookup(q)==0) break;
            if(c==1) goto a;
            if(a<=M) {
                if((a&1)==0) Q.install(a>>1);
                else if(a<=L) Q.install(a);
            }
        }
        if((a = GCD(q,b)) > 1) return a;
        for(a=q, c*=q;;) {
            r+=s;
            qr = ldiv(r,a);
            r=s-r;
            if(r==b) {
                if((a&1)==0) a>>=1;
                return a;
            }
            b-=r;
            b*=q;
            c+=b;
            b=r;
            r+=s;
            qr = ldiv(r,c);
            r=s-r;
            if(r==b) {
                if((c&1)==0) c>>=1;
                return c;
            }
            b-=r;
            b*=q;
            a+=b;
            b=r;
        }
    }
a:  ;
    PrimeSeq ps; ps.next();
    long p;
b:  ;
    if(bit(n,1))
        do p = ps.next(); while((p&3)!=3);
    else do p = ps.next(); while((p&3)!=1);
    if(n%p==0) return p;
    mul(D,n,p);
    SqrRoot(d,D);
    if(NumBits(d) > SQUFOF_MAX_BITS)
        Error("squfof: n is too large");
    conv(s,d);
    a=1; b=s; if((b&1)==0) b--;
    conv(d,b); d*=b; sub(d,D,d); d>>=2;
    conv(c,d);
    L = SqrRoot(s); M=L*p;
    Q.clear();
    s++;
    s>>=1; t=s-1;
    b>>=1; r=b;
    for(;;) {
        r+=s;
        qr = ldiv(r,c);
        r=t-r;
        b-=r;
        b*=q;
        a+=b;
        b=r;
        if(a==1) goto b;
        if(c<=M) {
            if(c%p==0) Q.install(c/p);
            else if(c<=L) Q.install(c);
        }
        r+=s;
        qr = ldiv(r,a);
        r=t-r;
        b-=r;
        b*=q;
        c+=b;
        b=r;
        if((q=IsSquare(a)) && Q.lookup(q)==0) break;
        if(c==1) goto b;
        if(a<=M) {
            if(a%p==0) Q.install(a/p);
            else if(a<=L) Q.install(a);
        }
    }
    if((a = GCD(q, (b<<1)|1)) > 1) return a;
    for(a=q, c*=q;;) {
        r+=s;
        qr = ldiv(r,a);
        r=t-r;
        if(r==b && a!=p) {
            if(a%p==0) a/=p;
            return a;
        }
        b-=r;
        b*=q;
        c+=b;
        b=r;
        r+=s;
        qr = ldiv(r,c);
        r=t-r;
        if(r==b && c!=p) {
            if(c%p==0) c/=p;
            return c;
        }
        b-=r;
        b*=q;
        a+=b;
        b=r;
    }
}
