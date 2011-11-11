function [flag n] = isSquare(A,throwError)
%ISSQUARE
% Checks whether a given matrix is square and returns the corresponding
% boolean value. If specified in the input, throws an error if the matrix
% is not square.
%
% Syntax:
%   [flag n] = isSquare(A,throwError)
%
% Input arguments:
% A             : Any matrix
% throwError    : Whether to throw an error in case the matrix is not
%                 square.
%
% Output arguments
% flag          : A boolean flag indicating whether the given matrix "A" is
%                 square or not.
% n             : The size of the matrix if the matrix is square or the
%                 number or rows otherwise.
%
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