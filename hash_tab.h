#include<NTL/vector.h>
using namespace NTL;
#define TAB_SIZE_MIN 1
#define TAB_SIZE_MAX 24

template<class T>
struct HashTable {
    long n,l,i,j;
    Vec<Vec<T> > elem;
    HashTable(long l1=0) : l(l1) {
        if(l < TAB_SIZE_MIN) l = TAB_SIZE_MIN;
        else if(l > TAB_SIZE_MAX) l = TAB_SIZE_MAX;
        n = (1<<l);
        elem.SetLength(n);
    }
    long hash_val(const T& t) {
        conv(j,t);
        for(i=0; j; j>>=l) i ^= j&(n-1);
        return	i;
    }
    void clear() {
        elem.kill();
        elem.SetLength(n);
    }
    T* lookup(const T& t) {
        i = hash_val(t);
        for(j=0; j < elem[i].length(); j++)
            if(elem[i][j] == t) break;
        if(j == elem[i].length()) return 0;
        return &(elem[i][j]);
    }
    inline void install(const T& t) {
        elem[hash_val(t)].append(t);
    }
};

