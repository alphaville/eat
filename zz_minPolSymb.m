function poly = zz_minPolSymb(A)
%zz_minPolySymb
%
% Uses MuPad to calculate the minimum polynomial symbolically
%
% Example...
% A=round(rand(5,5)*10); % A random 5-by-5 square matrix
% poly = zz_minPolySymb(A); % The minimal polynomial of A
% error=norm(polymatrixval(poly,A),'fro'); 
%
% Warning:
% The underlying algorithm is not numerically stable.

poly = sym2poly(evalin(symengine, ['map(poly2list(linalg::minpoly(' char(sym(A)) ...
    ',x)),c->c[1])']));