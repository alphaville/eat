function V = polymatrixval(poly,A)

polyLength=size(poly,2);
n=size(A,1);
assert(size(A,2)==n,'Input inconsistency: A is not a square matrix');

T=A;
V = poly(polyLength)*eye(n);
for i=1:polyLength-1
    V = V + poly(polyLength-i)*T;
    T = T*A;
end