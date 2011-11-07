function [polynomial distinctEvals exponents] = zz_minPolAlg(A,tol)
%zz_minPolAlg
%
% Returns the minimum polynomial of A using solely algebraic methods. This
% method is not numerically stable since it depends on the rank() of a
% matrix
%

% Step 1: Collect all distinct eigenvalues of A
n=size(A,1);
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
distinctEvals = evals(1);
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
        distinctEvals = [distinctEvals;toRegister];
        sizeDistinct = sizeDistinct+1;
    end
end