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
poolsize=matlabpool('size');
if poolsize==0
    disp('MATLAB pool was closed - now starting');
    matlabpool open;
end
parfor i=1:N
    VertexGNB{i}=polyOverAppr1(A,t1+(i-1)*h,t1+i*h);
end

Vertex={};
s=length(VertexGNB{1});
Vertex=[Vertex VertexGNB{1}{1}];
for i=1:N
    for j=1:s
        matrix = VertexGNB{i}{j};
        if (~contains(matrix,Vertex))
            Vertex = [Vertex matrix];
        end
    end
end

function flag = contains(matrix,Vcell)
n=length(Vcell);
flag=false;
if n==0
    flag=false;
else   
    for i=1:n
        if norm(matrix-Vcell{i})<1E-6
            flag=true;
            break;
        end
    end
end