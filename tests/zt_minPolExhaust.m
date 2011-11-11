N=5;
A=round(rand(N,N)*10*(-1)).*(-1).^round(rand(N,N)*2.+1);
poly1=zz_minPolExhaust(A);
errorExhaust=norm(polymatrixval(poly1,A));
assert(errorExhaust<=1E-6);
poly2=zz_minPolSymb(A);
errorSymb=norm(polymatrixval(poly2,A));
assert(errorSymb<=1E-6);