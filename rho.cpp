#include<NTL/ZZ.h>
using namespace NTL;

static long GCD_INTERVAL = 100;

void	pollard_rho(ZZ& d, const ZZ& n) {
	long	a(1),i;
a:	ZZ	u(2),v(2),q(1),s,t;
	do {	s=u; t=v;
		for(i=0; i<GCD_INTERVAL; i++) {
			SqrMod(u,u,n); AddMod(u,u,a,n);
			SqrMod(v,v,n); AddMod(v,v,a,n);
			SqrMod(v,v,n); AddMod(v,v,a,n);
			SubMod(d,u,v,n);
			MulMod(q,q,d,n);
		}
		GCD(d,q,n);
	} while(IsOne(d));
	if(d<n) return;
	do {	SqrMod(s,s,n); AddMod(s,s,a,n);
		SqrMod(t,t,n); AddMod(t,t,a,n);
		SqrMod(t,t,n); AddMod(t,t,a,n);
		sub(q,s,t);
		GCD(d,q,n);
	} while(IsOne(d));
	if(d<n) return;
	a++; goto a;
}

void	brent_rho(ZZ& d, const ZZ& n) {
	long	a(1),r,i,j;
a:	ZZ	u(2),q(1),s,t;
	for(r=1; r>0; r<<=1) {
		s=u;
		for(i=0; i<r; i++) {
			SqrMod(u,u,n);
			AddMod(u,u,a,n);
		}
		for(i=j=0; i<r;) {
			t=u;
			j += GCD_INTERVAL;
			if(j>r) j=r;
			for(; i<j; i++) {
				SqrMod(u,u,n);
				AddMod(u,u,a,n);
				SubMod(d,s,u,n);
				MulMod(q,q,d,n);
			}
			GCD(d,q,n);
			if(!IsOne(d)) goto b;
		}
	}
b:	if(d<n) return;
	do {	SqrMod(t,t,n);
		AddMod(t,t,a,n);
		sub(q,s,t);
		GCD(d,q,n);
	} while(IsOne(d));
	if(d<n) return;
	a++; goto a;
}
