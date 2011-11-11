function [Vertex Details] = hullJordIntegral(A,B,t1,t2)
%hullJordIntegral calculates a polytopic overapproximation for functions of
%the form g(r)=Integral [ expm(A*(r-s))B ds ] for r in [t1,t2] where
%usually t1=0 and t2>0. A is a squarre n-by-n matrix and B is a m-by-n one.
%Symbolic calculations are involved
%
% Syntax [Vertex Details] = hullJordIntegral(A,B,t1,t2)
%    Input:
%      A is a square matrix (n-by-n)
%      B is an n-by-m matrix. If B is the unit matrix (B=eye(n)) then the
%      result is the hull of g(r)=Integral [ expm(A*(r-s)) ds ].
%      t1,t2 define a time interval [t1,t2] which stands for the domain of
%      the function g(r)=Integral [ expm(A*(r-s))B ds ].
%    Output:
%      Vertex is a cell containing the various vertices of the polytopic
%      overapproximation of the range of f over [t1,t2].
%      Details is a cell of dimensions 1-by-K whose entries are structures.
%      Each structure has the following fields:
%        - fun: The function gamma_i as in equation (22) in [1]. That is
%        the function f(t)=expm(A*t) is written as follows:
%            f(t)=SUM_{i=1,...,K}gamma_i(t)*S_i
%        - S: A 0-1 matrix as above
%        - minFun and maxFun are the minimum and maximum values of the
%        function `fun` over the interval [t1,t2] calculated by means of
%        the MATLAB function exhaustiveSearch.
%
%   References:
%   [1] W.P.M.H. Heemels, N. vd Wouw, et al, "Comparison of
%   Overapproximation Methods for Stability Analysis of Networked Control
%   Systems"

syms r s;
n=size(A,1);
m=size(B,2);
% Symbolically calculate the integral
% g(r) = INT { expm(A*(r-s))*B ds } from 0 to r
% The dimensions of the result are m-by-n
symbolic = simple(int(expm(A*(r-s))*B,s,0,r));

% k is an index used to enumerate the number of non-zero function inside
% symbolicExpo.
k=0;

for i=1:n
    for j=1:m
        if (symbolic(i,j)~=0)
            % When k==0 no value has been registered to the `Details`
            % structure. In this case the function is registered and the
            % matrix S has only one 1 and all other entries are set to 0. 
            if (k==0)
                k=k+1;
                Details{k}.fun=symbolic(i,j);
                S = zeros(n,m);
                S(i,j)=1;
                Details{k}.S=S;
            end
            
            % Have we already registered this function in `Details` ?
            % Search throughout the entries of Details:
            if (k>0)
                isRegistered=0;
                for p=1:k
                    if (symbolic(i,j)==Details{p}.fun) %found already registered!
                        isRegistered=1;
                        S = Details{p}.S;
                        S(i,j)=1;
                        Details{p}.S=S;
                    end
                end
                % In the following case no function was found in `Details`
                % matching the candidate one, so it is registered.
                if (~isRegistered)
                    k=k+1;
                    Details{k}.fun=symbolic(i,j);
                    S = zeros(n,m);
                    S(i,j)=1;
                    Details{k}.S=S;
                end
            end
        end
    end
end

% Now the minimum and maximum values of Details{*}.fun are calculated and
% are stored in Details as .minFun and .maxFun. For this purpose the
% function exhaustiveSearch is used.
for i=1:k
    [mm MM]=exhaustiveSearch(Details{i}.fun,t1,t2,50);
    Details{i}.minFun=mm;
    Details{i}.maxFun=MM;
end

% Form the cartesian product space of min and max...
X=cell(k,1);
for i=1:k
    X{i}={Details{i}.minFun,Details{i}.maxFun};
end
PROD = cartprod(X);
nProd = size(PROD,2);

Vertex=cell(nProd,1);
for q=1:nProd
   vertice = zeros(n,m);
   for i=1:k
       tmp = PROD{q}(i);       
       vertice = vertice + tmp{1}*Details{i}.S; 
   end
          Vertex{q}=vertice;
end