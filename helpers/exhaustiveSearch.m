function [m M]=exhaustiveSearch(fun,t1,t2,N)
%EXHAUSTIVESEARCH 
%exhaustiveSearch is a function that returns the minimum and maximum values
%of a function of one variable over a specified closed interval [t1,t2] by
%means of an exhaustive search algorithm. This provides global extreme
%values.
%
% Syntax [min max]=exhaustiveSearch(fun,t1,t2,N)
%        [min max]=exhaustiveSearch(fun,t1,t2)
% Input:
%   fun is the function for which minimum and maximum values are calculated
%   t1 and t2 are the lower and upper bounds of the interval over which the
%   etreme values are calculated
%   N is the number of points used for the partition of the interval
%   [t1,t2]. If not specified, defaults to 1000.
if (nargin==3)
    N=1000;
end
h=(t2-t1)/N;
m=subs(fun,t1);
M=m;
for i=1:N
    temp = subs(fun,t1+i*h);
    M = max([M,temp]);
    m = min([m;temp]);
end