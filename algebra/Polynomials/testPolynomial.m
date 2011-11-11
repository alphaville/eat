%Test the zero polynomial
p=polynomial();
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

% Test Copying
p=randonomial(4);
q=randonomial(p);
assert(p==q);
q=3*q;
assert(p~=q);
p=randonomial(4);
q=p;
q=3*q;
assert(p~=q);

% Various tests
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

%Test random
rand1=randonomial();
rand2=randonomial(3);
assert(deg(rand2)==3);
plot(rand1);
hold on;
plot(rand2,'r');

%Test minimal polynomial
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

% Test minimal polynomial again
A=[ 1 1 1
    0 1 1
    0 0 1 ];
p=minimalPolynomial(A);
assert(p(1)==0);
assert(deg(p)==3);
assert(p==Polynomial([1 -3 3 -1]));
assert(norm(p(A))==0);



   


  