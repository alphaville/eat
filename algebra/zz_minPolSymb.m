function poly = zz_minPolSymb(A,mode)
%zz_minPolySymb
%
% Uses MuPad to calculate the minimum polynomial symbolically
%
% Syntax:
% poly = zz_minPolSymb(A,mode)
%
%
% Example...
% A=round(rand(5,5)*10); % A random 5-by-5 square matrix
% poly = zz_minPolySymb(A); % The minimal polynomial of A
% error=norm(polymatrixval(poly,A),'fro'); 
%
% Warning:
% The underlying algorithm is not numerically stable.

error(nargchk(1,2,nargin,'string'));
isSquare(A,true);

if nargin==1
    mode='Numerical';
end
if ~strcmp(mode,'Numerical') && ...
    ~strcmp(mode,'Rational') &&...
    ~strcmp(mode,'Real')  &&...
    ~strcmp(mode,'Natural')
    error('Wrong Mode: Choose between "Rational", "Real" and "Numerical"!');
end
symMatrix=['Dom::Matrix(Dom::' mode ')(' char(sym(A)) ')'];
symbolic=['poly2list(linalg::minpoly(' symMatrix ',x))']
polynomial = evalin(symengine, symbolic)
pMaxTerm=polynomial(1);
maxExpo=eval( pMaxTerm(2));
poly=zeros(1,maxExpo+1);
for i=1:length(polynomial)
    polyTerm=eval(polynomial(i));
    poly(maxExpo-polyTerm(2)+1)=polyTerm(1);
end