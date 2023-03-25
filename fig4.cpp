#include<NTL/ZZ.h>
#include<fstream>
using namespace NTL;

long squfof(const ZZ&);
void brent_rho(ZZ&, const ZZ&);
void mpqs(ZZ&, const ZZ&);
void morbri(ZZ&, const ZZ&);

double timer_squ(const ZZ& p) {
    double t;
    t = GetTime(); squfof(p);
    t = GetTime() - t;
    return t;
}

double timer(void(*f)(ZZ&, const ZZ&), const ZZ& p) {
    double t;
    ZZ q;
    t = GetTime(); f(q,p);
    t = GetTime() - t;
    return t;
}

enum { squ_idx, rho_idx, mpqs_idx, mb_idx };

void fig4(long idx, const char *fname, long N, long l1, long l2, long L=0) {
    std::ofstream ofs(fname);
    long i,l;
    double t,t1;
    ZZ p,q;
    for(l=l1; l<=l2; l++) {
        for(i=0, t=0; i<N; i++) {
            if(l<L) GenPrime(p,L-l);
            else GenPrime(p,l);
            do GenPrime(q,l); while(p==q);
            p*=q;
            switch(idx) {
                case squ_idx:  t += timer_squ(p);        break;
                case rho_idx:  t += timer(brent_rho, p); break;
                case mpqs_idx: t += timer(mpqs, p);      break;
                case mb_idx:   t += timer(morbri, p);    break;
            }
        }
        t /= N;
        std::cout << l << '\t' << t << std::endl;
        ofs << l << '\t' << t << std::endl;
    }
}

main() {
    fig4(squ_idx,  "fig4a.txt", 5,  5, 60);
    fig4(rho_idx,  "fig4b.txt", 5,  5, 55);
    fig4(mpqs_idx, "fig4c.txt", 10, 20, 60);
    fig4(mb_idx,   "fig4d.txt", 10, 20, 60);
    fig4(squ_idx,  "fig5a.txt", 20, 3,  57, 60);
    fig4(rho_idx,  "fig5b.txt", 20, 3,  57, 60);
}
