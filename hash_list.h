#include<list>
#define	HLIST_SIZ_MIN   1
#define	HLIST_SIZ_MAX   24
#define	HLIST_DEPTH_MAX 64

template<class T>
struct HashList {
    long n,l,d,i,j;
    std::list<T> *elem;
    typename std::list<T>::iterator p;
    HashList(long l1=0, long d1=0) : l(l1), d(d1) {
        if(l < HLIST_SIZ_MIN) l = HLIST_SIZ_MIN;
        else if(l > HLIST_SIZ_MAX) l = HLIST_SIZ_MAX;
        if(d<=0 || d > HLIST_DEPTH_MAX) d = HLIST_DEPTH_MAX;
        n = (1<<l);
        elem = new std::list<T>[n];
    }
    ~HashList() { delete[] elem; }
    long hash_val(const T& t) {
        conv(j,t);
        for(i=0; j; j>>=l) i ^= j&(n-1);
        return	i;
    }
    inline void clear() {
        for(i=0; i<n; i++) elem[i].clear();
    }
    T* lookup_or_install(const T& t) {
        i = hash_val(t);
        for(p = elem[i].begin();
            p!= elem[i].end(); p++)
            if(!(t>*p)) break;
        if(t==*p) return &(*p);
        else if(elem[i].size()<d) elem[i].insert(p,t);
        return 0;
    }
};

