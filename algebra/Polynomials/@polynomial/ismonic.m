function t = ismonic(obj)
%ISMONIC
% Checks whether the polynomial is monic, i.e. whether the high-order
% coefficient of the polynomial is approximately 1. Approximately is meant
% up to a tolerance of 1E-9.
%
% Syntax:
% t = ismonic(obj)
%
% Input arguments:
%  obj : A polynomial object
% 
% Output arguments:
%  t   : A boolean flag which is true if and only if the input polynomial
%        is monic.
error(nargchk(1,1,nargin));
assert(ispolynomial(obj));
if deg(obj)==-1
    t=0;
    return;
end
coeffs=obj.coeff;
t=abs(coeffs(1)-1)<1E-9;
