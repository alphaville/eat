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

error(nargchk(1,1,nargin,'This function requires 1 input argument'));

isSquare(A,true);
symbolic=['poly2list(linalg::minpoly(' char(sym(A)) ',x))'];
polynomial = evalin(symengine, symbolic);
pMaxTerm=polynomial(1);
maxExpo=eval( pMaxTerm(2));
poly=zeros(1,maxExpo+1);
for i=1:length(polynomial)
    polyTerm=polynomial(i);       
    poly(maxExpo-eval(polyTerm(2))+1)=eval(polyTerm(1));
end