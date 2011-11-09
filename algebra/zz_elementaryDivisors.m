function [V nMulti]= zz_elementaryDivisors(A,tol)

[~,n]=isSquare(A,true);
if nargin==1
    tol=1E-7;
end

[distinctEvals sizeDistinct] = zz_distinctEvals(A,tol);
V=[];

% Clear conjugate eigenvalues from distinctEvals
distinctEvalsClean = zeros(sizeDistinct,2);
distinctEvalsClean(1,:)=distinctEvals(1,:);
cleanSize=1;
for i=2:sizeDistinct
    isAlreadyRegistered = 0;
    currentEval=distinctEvals(i,1);    
    for j=1:i-1
        if abs(currentEval-conj(distinctEvalsClean(j)))<=tol ||...
                abs(currentEval-distinctEvalsClean(j))<=tol
            isAlreadyRegistered=1;
        end
    end
    if ~isAlreadyRegistered
        cleanSize = cleanSize + 1;
        distinctEvalsClean(cleanSize,1) = currentEval;
        distinctEvalsClean(cleanSize,2) = distinctEvals(i,2);
    elseif isreal(currentEval)
        disp('!!!!');
        distinctEvalsClean(cleanSize,2) = distinctEvalsClean(cleanSize,2)+1;
    end
end

% Create elementary divisors
nMulti=1;
for i=1:cleanSize
    currentEval=distinctEvalsClean(i,1);
    multi=distinctEvalsClean(i,2);
    nMulti=nMulti*multi;
    if isreal(currentEval)
        V=[V struct('polynomial',[1 -currentEval],'multiplicity',multi)];
    else
        re=real(currentEval);
        im=imag(currentEval);
        V=[V struct('polynomial',[1 -2*re re^2+im^2],'multiplicity',multi)];
    end
end