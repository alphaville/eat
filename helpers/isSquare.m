function [flag n] = isSquare(A,throwError)
if nargin==1
    throwError=false;
end
n=size(A,1);
m=size(A,2);
flag=n==m;
if throwError
    assert(flag,...
        ['Matrix not square:: Dimensions :',num2str(n),'-by-',num2str(m)]);
end