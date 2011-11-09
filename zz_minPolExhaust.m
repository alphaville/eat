function [poly Rcombinations] = zz_minPolExhaust(A,tol)

%assign default tolerance
if nargin==1
    tol=1E-7;
end
isSquare(A,true);%check whether A is square - throw error if not
% Calculate elementary divisors
% of the characteristic polynomial of A
[elDiv nMulti]= zz_elementaryDivisors(A,tol);

if nMulti==1
    % There are only single eigenvalues
    poly=zz_charPolyAlg(A,tol);
else
    % There are multiple eigenvalues...
    vLength=length(elDiv);
    multiCell=cell(1,vLength);
    for i=1:vLength
        m=elDiv(i).multiplicity;
        mCell=cell(1,m);
        for j=1:m
            mCell{j}=j;
        end
        multiCell{i}=mCell;
    end
    Rcombinations=cartprod(multiCell); %possible combinations of exponents
    minDim=Inf;
    poly=-1;
    for i=1:length(Rcombinations)
        testPoly=candidatePoly(elDiv,Rcombinations{i});
        criterion=norm(polymatrixval(testPoly,A),'fro');
        dim=length(testPoly);
        if dim<minDim && criterion<tol
            minDim=dim;
            poly=testPoly;
        end
    end
end
if poly==-1
    disp('No polynomial found - Not even the char. pol. does annihilate the matrix!');
    disp('This is due to numerical errors - Try again with different tolerance');
end

function p=candidatePoly(elDivisors,comb)
%CANDIDATEPOLY
%
% Given the elementary divisors v1,v2,...vs, create a polynomial
% v1^c1*v2^c2*...vs^cs using the exponents in comb
p=1;
if iscell(comb)
    for i=1:length(elDivisors)
        p=conv(p,expose(elDivisors(i).polynomial,comb{i}));
    end
else
    p=conv(p,expose(elDivisors(1).polynomial,comb));
end

function pol=expose(polynomial, exponent)
%EXPOSE
%
% Expose a polynomial to a given power.
pol=1;
for i=1:exponent
    pol=conv(pol,polynomial);
end