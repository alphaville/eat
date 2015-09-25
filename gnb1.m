function Vertex = gnb1(A,t1,t2,N)
%GNB1
% calculates a polytopic overapproximation for functions of
% the form f(t)=expm(A*t) for t in [t1,t2] for given t1 < t2 using the
% gridding and bounding algorithm with an underlying Jordan-based
% overapproximation as described in [1].
%
% Syntax Vertex = gnb1(A,t1,t2,N)
%    Input:
%      A is a square matrix (n-by-n)
%      t1,t2 define a time interval [t1,t2] which stands for the domain of
%      the function f(t)=expm(A*t).
%      N is the number of subintervals used by the Gridding and Bounding
%      approximation technique. Defaults to 100 unless otherwise
%      specified.
%    Output:
%      Vertex is a cell containing the various vertices of the polytopic
%      overapproximation of the range of f over [t1,t2].
%
%   References:
%   [1] W.P.M.H. Heemels, N. vd Wouw, et al, "Comparison of
%   Overapproximation Methods for Stability Analysis of Networked Control
%   Systems"
if (nargin==3)
    N=100;
end
h=(t2-t1)/N;
VertexGNB=cell(N,1);
% start MATLAB parallel processing pool:
parfor i=1:N
    VertexGNB{i}=polyOverAppr1(A,t1+(i-1)*h,t1+i*h);
end

s=length(VertexGNB{1});
Vertex={  real(double(VertexGNB{1}{1}))   };

for i=1:N
    for j=1:s
        matrix = VertexGNB{i}{j};
        if (~z_contains(matrix,Vertex))
            Vertex{end+1} = real(double(matrix));
        end
    end
end

function flag = z_contains(matrix,Vcell)
n=length(Vcell);
flag=false;
if n==0
    flag=false;
else
    for i=1:n
        if norm(double(matrix)-double(Vcell{i})) < 1E-6
            flag=true;
            break;
        end
    end
end