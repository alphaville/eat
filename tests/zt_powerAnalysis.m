clear all;
disp('Clear - GO!');

X=[ -1   4
    -2  -1 ];
Y=[ 1   1
    0   1 ];
A=blkdiag(X,Y);
q=minimalPolynomial(A);
v=deg(q); % dimension of the minimal polynomial
for i=1:3*v
    p=zz_powerAnalysis(A,i,q);
    assert(norm(p(A)-A^i)<1E-6);
end
disp('Test 1 :: passed [OK]');

clear q v p ;
A=blkdiag(Y,0,Y,ones(4,4),0);
q=minimalPolynomial(A);
v=deg(q); % dimension of the minimal polynomial
for i=1:3*v
    p=zz_powerAnalysis(A,i,q);
    assert(norm(p(A)-A^i)<1E-6);
end
disp('Test 2 :: passed [OK]');

A=[1 2 4;3 4 -1;0 1 0];
q=minimalPolynomial(A);
v=deg(q); % dimension of the minimal polynomial
for i=1:3*v
    p=zz_powerAnalysis(A,i,q);
    assert(norm(p(A)-A^i)<1E-6);
end
disp('Test 3 :: passed [OK]');

A=[1 -3 -5
    100 2 1
    1 0 -23];
A=A/norm(A,'fro');
q=minimalPolynomial(A);
v=deg(q); % dimension of the minimal polynomial
for i=1:3*v
    p=zz_powerAnalysis(A,i,q);
    assert(norm(p(A)-A^i)<1E-7);
end
disp('Test 4 :: passed [OK]');

A=blkdiag(X,Y,Y,Y,Y,X,X);
q=minimalPolynomial(A);
v=deg(q); % dimension of the minimal polynomial
for i=1:3*v
    p=zz_powerAnalysis(A,i,q);
    assert(norm(p(A)-A^i)<1E-6);
end
disp('Test 5 :: passed [OK]');