function [distinctEvals sizeDistinct] = zz_distinctEvals(A,tol)
%zz_distinctEvals
%
%
%
%

if nargin==1
    tol=1E-7;
end
[~,n]=isSquare(A,true);
evals=eig(A);
distinctEvals=ones(n,2);
distinctEvals(1,1)=prv_shapeUpComplex(evals(1),tol);

% iterate over all eigenvalues of A
sizeDistinct=1;
for i=2:n
    isAlreadyRegistered=0;
    currentEval=prv_shapeUpComplex(evals(i),tol);
    
    % Iterate over already registered eigenvalues in "distinctEvals"
    for j=1:sizeDistinct
        if abs(currentEval-distinctEvals(j))<=tol
            isAlreadyRegistered=1;
            distinctEvals(j,2)=distinctEvals(j,2)+1;
            break;
        end
    end
    % Register new eigenvalue
    if isAlreadyRegistered==0
        sizeDistinct = sizeDistinct +1;
        distinctEvals(sizeDistinct,1)=currentEval;
    end
end
distinctEvals=distinctEvals(1:sizeDistinct,:);

function newNumber = prv_shapeUpComplex(oldNumber,tol)
% Treat false complex numbers:
newNumber=oldNumber;
if abs(imag(oldNumber))<=tol
    newNumber=real(oldNumber);
end
if abs(real(oldNumber))<=tol
    newNumber=imag(oldNumber)*1i;
end