classdef RandomPolynomial < Polynomial
    methods
        function obj = RandomPolynomial(varargin)
            if nargin==0
                deg=round(rand()*10);
            elseif nargin>=1
                if isa(varargin{1},'Polynomial')
                    %Copy another polynomial
                    obj.coeff = varargin{1}.coeff;
                    return;
                end
                deg=varargin{1};
            end
            coefficients=rand(1,deg+1);
            if nargin==2 && strcmp(varargin{2},'integer')
                coefficients=round(coefficients*10);
            end
            obj.coeff=coefficients;
            
        end
    end
end