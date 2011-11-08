function [k poly error] = zz_minPolGramSchmidt(A,tol)
%zz_minPolGramSchmidt
%
% Uses the Gram-Schmidt algorithm to calculate the minimal polynomial of a
% given matrix. 
%
% WARNING: This algorithm, while in exact arithmetic performs well, it
% fails to give any results in most real cases. The algorithm is actually
% so ill-conditioned that you might even have difficulty finding a single
% working example. For that, consider:
% 
% A=diag([1;1;2;2;3]);
% [k poly error] = zz_minPolGramSchmidt(A,tol);
%
% The error value returned is 1.0993e-12. However, if you replace A with
% rand(4,4), it will not work.
%
% Syntax:
% [k U]=zz_minPolGramSchmidt(A[, tol])
%
% Input...
% A     :   An n-square matrix
% tol   :   Tolerance. Used to tell whether two eigenvalues are equal
%
% Output...
% k     :   The degree of the minimal polynomial of A
% poly  :   The minimal polynomial of A
% error :   The norm of poly(A) (should be 0)
%
%See also:
% polymatrixval, zz_minPolAlg


% Check input
if nargin==1
    tol=1E-7;
end
n=size(A,1);
assert(size(A,2)==n,'Input inconsistency: A is not a square matrix');

% Gram-Schmidt algorithm
U=cell(1);
v=zeros(1,n);
normFroEye=sqrt(n);
U{1}=eye(n)/normFroEye;
v(1)=normFroEye;
AJ=A;
k=1;
for j=1:n
    %... assert(norm(AJ-A^j,'fro')<1E-8);
    Sum=zeros(n);
    for i=1:j
        Sum = Sum + trace(U{i}*AJ)*U{i};
    end
    nextMatr=AJ-Sum;
    normU=norm(nextMatr,'fro');
    if (normU <= tol)
        break;
    else
        U{i+1}=nextMatr/normU;
        v(j+1)=normU;
        k=k+1;
    end
    AJ=AJ*A;
end

v=v(1:k);
r=zeros(k,k);
for i=1:k
    for j=i:k+1
        r(i,j)=trace(U{i}*A^(j-1));
    end
end


% Assertions
assert(k<=n+1);
assert(norm(v(1)*U{1}-eye(n))<1E-7);
for i=1:k-1
    Sum=v(i+1)*U{i+1};
    for j=1:i
        Sum = Sum + r(j,i+1)*U{j};
    end
    assert(norm(A^i-Sum)<1E-7);
end

c=r(:,end);
R=zeros(k,k);
for i=1:k
    for j=i:k
        if i==j
            R(i,i)=v(i);
        else
            R(i,j)=r(i,j);
        end
    end
end
poly=[1 fliplr(-(R\c)')];
error = norm(polymatrixval(poly,A),'fro');