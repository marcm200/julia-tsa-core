-------------------
Cubic: z := z^3 + c
-------------------

inline void getBoundingBoxfA(PlaneRect& A,PlaneRect& fA) {
	fA.x0=A.x0*A.x0*A.x0-(3*maximumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedCre;
	fA.x1=A.x1*A.x1*A.x1-(3*minimumD(A.x0*minimumD(A.y0*A.y0,A.y1*A.y1),A.x0*maximumD(A.y0*A.y0,A.y1*A.y1),A.x1*minimumD(A.y0*A.y0,A.y1*A.y1),A.x1*maximumD(A.y0*A.y0,A.y1*A.y1)))+seedCre;
	fA.y0=3*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y1*A.y1*A.y1)+seedCim;
	fA.y1=3*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0,A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0,A.x1*A.x1)*A.y1)-(A.y0*A.y0*A.y0)+seedCim;
}

--------------------------------------
Quartic: Bounding box for z := z^4 + c
--------------------------------------

inline void getBoundingBoxfA(PlaneRect& A,PlaneRect& fA) {
	fA.x0=minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedCre;
	fA.x1=maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)-(6*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),minimumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1)))+maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1)+seedCre;
	fA.y0=4*minimumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*maximumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedCim;
	fA.y1=4*maximumD((A.x0*A.x0*A.x0)*A.y0,(A.x0*A.x0*A.x0)*A.y1,(A.x1*A.x1*A.x1)*A.y0,(A.x1*A.x1*A.x1)*A.y1)-(4*minimumD(A.x0*(A.y0*A.y0*A.y0),A.x0*(A.y1*A.y1*A.y1),A.x1*(A.y0*A.y0*A.y0),A.x1*(A.y1*A.y1*A.y1)))+seedCim;
}

--------------------
Pentic: z := z^5 + c
--------------------

inline void getBoundingBoxfA(PlaneRect& A,PlaneRect& fA) {
	fA.x0=A.x0*A.x0*A.x0*A.x0*A.x0-(2*(5*maximumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*minimumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedCre;
	fA.x1=A.x1*A.x1*A.x1*A.x1*A.x1-(2*(5*minimumD((A.x0*A.x0*A.x0)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x0*A.x0*A.x0)*maximumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*minimumD(A.y0*A.y0,A.y1*A.y1),(A.x1*A.x1*A.x1)*maximumD(A.y0*A.y0,A.y1*A.y1))))+5*maximumD(A.x0*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x0*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*minimumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1),A.x1*maximumD(A.y0*A.y0*A.y0*A.y0,A.y1*A.y1*A.y1*A.y1))+seedCre;
	fA.y0=5*minimumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*maximumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y0*A.y0*A.y0*A.y0*A.y0+seedCim;
	fA.y1=5*maximumD(minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,minimumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y0,maximumD(A.x0*A.x0*A.x0*A.x0,A.x1*A.x1*A.x1*A.x1)*A.y1)-(2*(5*minimumD(minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),minimumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y0*A.y0*A.y0),maximumD(A.x0*A.x0,A.x1*A.x1)*(A.y1*A.y1*A.y1))))+A.y1*A.y1*A.y1*A.y1*A.y1+seedCim;
}

