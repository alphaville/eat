function [polynomial distinctEvals exponents] = zz_minPolAlg(A,tol)
%zz_minPolAlg
%
% Returns the minimum polynomial of A using solely algebraic methods. This
% method is not numerically stable since it depends on the rank() of a
% matrix. To be used only for educational purposes and with small matrices.
% Afterwards you may use polymatrixval to check your result. Run:
%
% %Example:
% Z=.... % Define you matrix
% poly=zz_minPolAlg(Z);
% error=norm(polymatrixval(poly,A),'fro');
%
% Input:
% A     :   An n-square matrix
% tol   :   Tolerance. Used to tell whether two eigenvalues are equal

%First check input:
if nargin==1
    tol=1E-7;
end
n=size(A,1);
assert(size(A,2)==n,'Input inconsistency: A is not a square matrix');


% Step 1: Collect all distinct eigenvalues of A
[distinctEvals Ndistinct]= prv_distinctEvals(A,n,tol);

% Step 2: Determine the exponents
exponents=zeros(Ndistinct,1);
for i=1:Ndistinct
    exponents(i) = zz_eigenIndex(A,distinctEvals(i));
end


polynomial=1;
for i=1:Ndistinct
    if isreal(distinctEvals(i))
        polynomial = conv(polynomial,[1 -distinctEvals(i)]);
        for j=1:exponents(i)-1
            polynomial = conv(polynomial,[1 -distinctEvals(i)]);
        end
    else
        re=real(distinctEvals(i));
        img=imag(distinctEvals(i));
        p=[1 -2*re re^2+img^2];
        for j=1:exponents(i)
            polynomial = conv(polynomial,p);
        end
    end
    
end

function [distinctEvals sizeDistinct]= prv_distinctEvals(A,n,tol)

evals = eig(A);
distinctEvals = zeros(n,1);
distinctEvals(1) = evals(1);
sizeDistinct=1;
for i=1:n
    isAlreadyRegistered = 0;
    for j=1:sizeDistinct,
        if abs(evals(i)-distinctEvals(j))<=tol || abs(conj(evals(i))-distinctEvals(j))<=tol
            isAlreadyRegistered = 1;
        end
    end
    if ~isAlreadyRegistered
        toRegister=evals(i);
        if (imag(toRegister))<1E-7
            toRegister=real(toRegister);
        end
        distinctEvals(sizeDistinct+1) = toRegister;
        sizeDistinct = sizeDistinct+1;
    end
end
distinctEvals=distinctEvals(1:sizeDistinct);