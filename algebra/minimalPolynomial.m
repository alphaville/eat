function [poly verify err]= minimalPolynomial(A,options)
%MINIMALPOLYNOMIAL
%
% Calculate the minimal polynomial of a given matrix A by means of some
% subroutine of your choise. Available algorithms are:
%
% exhaust   :   Exhaustive search among the elementary divisors of the
%               characteristic polynomial of A.
% algebraic :   Algebraic method: The minimal polynomial is made out of the
%               eigenvalues of the matrix and their indices (uses the 
%               subroutine zz_eigenIndex.
% symbolic  :   Makes use of the symbolic algorithm for the calculation of
%               the minimal polynomial provided my MuPad.
% gs        :   Derivation of the minimal polynomial based on the
%               Gram-Schmidt orthogonalization algorithm.
%
% Syntax:
% [poly verify err]= minimalPolynomial(A,options)
%
% Input variables:
% A         :   An n-square matrix (mandatory).
% options   :   This is an optional input argument containing options for
%               the underlying method. It is a structure with the following
%               fields:
%   * tolerance : If the underlying routing uses a tolerance for the
%                 discrimination between eigenvalues of as a stop criterion
%                 this value will be used. If not provided, defaults to
%                 1E-7.
%   * algorithm : One of the available algorithms (exhaust, algebraic,
%                 symbolic, gs). If not specified, defaults to 'exhaust').
%   * verify    : Whether to verify that the calculated polynomial is
%                 indeed an annihilating polynomial for the given matrix.
%                 The norm ||p(A)|| is calculated and is exported to the
%                 output via the output argument err. The norm used is the
%                 Frobenius norm. If the error value is less than or equal
%                 to the specified tolerance, then the calculated
%                 polynomial annihilates A with respect to the specified
%                 tolerance.
%
% Output variables:

if nargin==1
    options=struct('tolerance',1E-7,'algorithm','exhaust','verify',false);
elseif nargin==2
    if ~isfield(options,'tolerance')
        options.tolerance=1E-7;
    end
    if ~isfield(options,'algorithm')
        options.algorithm='exhaust';
    end
    if ~isfield(options,'verify')
        options.verify=false;
    end
end

switch options.algorithm
    case 'exhaust'
        poly=zz_minPolExhaust(A,options.tolerance);
    case 'algebraic'
        poly=zz_minPolAlg(A,options.tolerance);
    case 'symbolic'
        poly=zz_minPolSymb(A);
    case 'gr'
        poly=zz_minPolGramSchmidt(A,tol);
    otherwise
        error(['Option : ',options.algorithm,' is not eligible for field "algorithm".']);
end
if options.verify
    
end
poly=Polynomial(poly);
err=norm(polymatrixval(poly,A),'fro');
verify=err<=options.tolerance;