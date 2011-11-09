function poly = zz_minPolExhaust(A)

[~,n]=isSquare(A,true);
V = zz_elementaryDivisors(A,tol);
poly=2*A;
