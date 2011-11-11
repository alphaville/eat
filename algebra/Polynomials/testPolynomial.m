p=Polynomial();
assert(p.coeff==0);
assert(p+p==p);
assert(p-p==p);
assert(+p==p);
assert(-p==p);
assert(100*p==p);
assert(p(1)==0);
A=[1 2;3 4];
assert(norm(p(A)-zeros(2,2))==0);
assert(deg(p)==-1);

p=Polynomial([1 2 3 4]);
assert( norm(p.coeff-[1 2 3 4])==0 );
assert( norm(-p.coeff+[1 2 3 4])==0 );
assert(p==p);
assert(p~=-p);
assert(strcmp('x^3 + 2*x^2 + 3*x + 4',char(p)));
assert(deg(p)==3);

p=Polynomial([1 2 3]);
q=Polynomial(p);
assert(p==q);

rand=RandomPolynomial();
rand=RandomPolynomial(10);
assert(deg(rand)==10);
plot(rand);


A =[
    4     2     0     0     0     0     0     0     0     0
    -2     4     0     0     0     0     0     0     0     0
    0     0     1     1     0     0     0     0     0     0
    0     0     0     1     1     0     1     0     0     0
    0     0     0     0     1     0     0     1     0     0
    0     0     0     0     0     4     2     0     0     0
    0     0     0     0     0    -2     4     0     0     0
    0     0     0     0     0     0     0     4     2     0
    0     0     0     0     0     0     0    -2     4     0
    0     0     0     0     0     0     0     0     0     1];
p=Polynomial(minimalPolynomial(A));
assert(norm(p(A))==0);
