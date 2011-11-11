classdef Polynomial
    %Polynomial Class
    %
    properties
        % coefficients of the polynomial in descending order. For example
        % the coefficients of the polynomial p=x^2+2x would be [1 2 0].
        % Trailing zeros in the beginning are automatically removed.
        coeff = 0;
    end
    
    methods
        function obj = Polynomial(args)
            % Constructor for a polynomial object with given coefficients.
            
            if nargin==0
                %No arguments => Zero polynomial
                obj.coeff=0;
                return;
            end
            if isa(args,'Polynomial')
                %Copy another polynomial
                obj.coeff = args.coeff;
                return;
            end
            if (isnumeric(args))
                obj.coeff=args;
            end
        end
        
        function flag=eq(obj1,obj2)
            %Whether two polynomials are equal (up to tolerance 1E-9)
            flag=norm(obj1.coeff-obj2.coeff,2)<=1E-9;
        end
        
        function flag=ne(obj1,obj2)
            %Whether two polynomials are equal (up to tolerance 1E-9)
            flag=norm(obj1.coeff-obj2.coeff,2)>1E-9;
        end
        
        function c = double(obj)
            % Converts the polynomial object to double.
            c = obj.coeff;
        end
        
        function str = char(obj)
            % Creates a formated display of the polynomial
            % as powers of x
            if all(obj.coeff == 0)
                str = '0';
                return;
            else
                d = length(obj.coeff)-1;
                s = cell(1,d);
                ind = 1;
                for a = obj.coeff;
                    if a ~= 0;
                        if ind ~= 1
                            if a > 0
                                s(ind) = {' + '};
                                ind = ind + 1;
                            else
                                s(ind) = {' - '};
                                a = -a;
                                ind = ind + 1;
                            end
                        end
                        if a ~= 1 || d == 0
                            if a == -1
                                s(ind) = {'-'};
                                ind = ind + 1;
                            else
                                s(ind) = {num2str(a)};
                                ind = ind + 1;
                                if d > 0
                                    s(ind) = {'*'};
                                    ind = ind + 1;
                                end
                            end
                        end
                        if d >= 2
                            s(ind) = {['x^' int2str(d)]};
                            ind = ind + 1;
                        elseif d == 1
                            s(ind) = {'x'};
                            ind = ind + 1;
                        end
                    end
                    d = d - 1;
                end
            end
            str = [s{:}];
        end
        
        function obj = set.coeff(obj,val)
            if ~isa(val,'double')
                error('Coefficients must be of class double')
            end
            ind = find(val(:).'~=0);
            if ~isempty(ind);
                obj.coeff = val(ind(1):end);
            else
                obj.coeff = val;
            end
        end
        
        function display(obj)
            % DISP Display object in MATLAB syntax
            c = char(obj); % char returns a cell array
            if iscell(c)
                disp(['     ' c{:}])
            else
                disp(c) % all coefficients are zero
            end
        end
        
        function b = subsref(a,s)
            % If p is a Polynomial and t is a numeric then p(t) is the
            % value of p at t. If A is a square n-by-n matrix (n>=2) then
            % p(A) is the evaluation of p at A. If B is a rectangular
            % matrix (not square, n-by-m) then an error message is thrown.
            if iscell(s.subs)
                A=s.subs{1};
                [flag n]=isSquare(A);
                if flag && n>=2
                    b=polymatrixval(a.coeff,A);
                    return;
                end
                
                if ~flag
                    error('p(B) is valid when B is a square matrix (or a number)');
                end
            end
            switch s(1).type
                case '()'
                    ind = s.subs{:};
                    b = a.polyval(ind);
                case '.'
                    switch s(1).subs
                        case 'coeff'
                            b = a.coeff;
                        case 'plot'
                            a.plot;
                        otherwise
                            if length(s)>1
                                b = a.(s(1).subs)(s(2).subs{:});
                            else
                                b = a.(s.subs);
                            end
                    end
                otherwise
                    error('Specify value for x as obj(x)')
            end
        end
        
        function r = plus(obj1,obj2)
            % Plus Implement obj1 + obj2 for Polynomials
            obj1 = Polynomial(obj1);
            obj2 = Polynomial(obj2);
            k = length(obj2.coeff) - length(obj1.coeff);
            r = Polynomial([zeros(1,k) obj1.coeff]+[zeros(1,-k) obj2.coeff]);
        end
        
        function r = minus(obj1,obj2)
            % MINUS Implement obj1 - obj2 for Polynomials
            obj1 = Polynomial(obj1);
            obj2 = Polynomial(obj2);
            k = length(obj2.coeff) - length(obj1.coeff);
            r = Polynomial([zeros(1,k) obj1.coeff]-[zeros(1,-k) obj2.coeff]);
        end
        
        function r = mtimes(obj1,obj2)
            % MTIMES Implement obj1 * obj2 for Polynomials
            obj1 = Polynomial(obj1);
            obj2 = Polynomial(obj2);
            r = Polynomial(conv(obj1.coeff,obj2.coeff));
        end
        
        function q = diff(obj,n)
            if nargin==1
                n=1;
            end
            % diff(obj) is the derivative of the DocPolynom obj
            c = obj.coeff;
            for i=1:n
                d = length(c) - 1;  % degree
                qc=c(1:d).*(d:-1:1);
                c = qc;
                if isempty(c)
                    break;
                end
            end
            q=Polynomial(qc);
            
        end
        
        function r = roots(obj)
            % roots(obj) returns a vector containing the roots of obj
            r = roots(obj.coeff);
        end
        
        function y = polyval(obj,x)
            % polyval(obj,x) evaluates obj at the points x
            y = polyval(obj.coeff,x);
        end
        
        function plot(obj)
            % plot(obj) plots the DocPolynom obj
            r = max(abs(roots(obj)));
            x = (-1.1:0.01:1.1)*r;
            y = polyval(obj,x);
            plot(x,y);
            title(['y = ' char(obj)])
            xlabel('X')
            ylabel('Y','Rotation',0)
            grid on
        end
        
        function y = uminus(obj)
            y = (-1)*obj;
        end
        
        function y = uplus(obj)
            y = obj;
        end
        
        function [y r] = mrdivide(obj1,obj2)
            [y r]=deconv(obj1.coeff,obj2.coeff);
        end
        
        function y = mpower(obj,x)
            y=Polynomial(1);
            for i=1:x
                y=y*obj;
            end
        end
        
        function d= deg(obj)
            if obj.coeff==0
               d=-1;
               return;
            end
            d = length(obj.coeff)-1;
        end
        
    end %end methods
end %end class
