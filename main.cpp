#include<NTL/ZZ.h>
using namespace NTL;

long squfof(const ZZ&);
void morbri(ZZ&, const ZZ&);

main() {
    ZZ p,q,n;
    long d;
    GenPrime(p,40);
    GenPrime(q,40);
    mul(n,p,q);
    std::cout << p << std::endl;
    std::cout << q << std::endl;
    std::cout << n << std::endl;
    d = squfof(n);
    std::cout << d << std::endl;
    morbri(p,n);
    std::cout << p << std::endl;
}
