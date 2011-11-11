function f = zz_powerAnalysis(A,z,q)
%zz_powerAnalysis
%
% Let A be a matrix and q its characteristic polynomial. Let v be the
% degree of the characteristic polynomial. Then A^j for all j>=0 can be
% written as a linear combination of I,A,A^2,...,A^(v-1). The result from
% this function is an object of class 'polynomial'.
%
% Syntax:
%  p = zz_powerAnalysis(A,z)
%  p = zz_powerAnalysis(A,z,q)
%
% Input Arguments:
% A : A square matrix
% z : An integer exponents
% q : The minimal polynomial of A (Optional). If not provided, it will be
%     calculated within this function.
%
% Output Arguments:
% p : A polynomial object so that A^z=p(A) and p has degree equal to the
%     degree of the minimal polynomial of A.

isSquare(A,true);
if nargin==2
    q=minimalPolynomial(A);
elseif nargin==3
    if isnumeric(q)
        q=polynomial(q);
    end
end

v=deg(q);

pol=-q.subpoly(2:v+1); % pol(A)=A^v
if z<=v-1
    f=polynomial([1 zeros(1,z)]);
    return;
elseif z==v
    f=pol;
    return;
end

xpol=pol;
for i=1:z-v
    mc=xpol{1};
    if deg(xpol)<v-1
        f=polynomial([1 0])*xpol;
        xpol=f;
    else
        s1=xpol.subpoly(2:v);
        s=polynomial([s1.coeff 0]);
        f=mc*pol+s;
        xpol=f;
    end
    
end

