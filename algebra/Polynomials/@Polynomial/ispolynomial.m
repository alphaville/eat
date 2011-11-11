function t = ispolynomial(obj)
error(nargchk(1,1,nargin));
t = strcmp('Polynomial',class(obj));