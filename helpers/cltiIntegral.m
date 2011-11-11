function [intVal Npoints] = cltiIntegral(x1,x2,x3,x4,x5)
%cltiIntegral
% Numerically calculates the integral of f(t)=expm(A*t)*B with respect to t
% over an interval [a,b]. The underlying numerical procedure is based on
% the trapezoid integration technique.
%
% Usage #1:
% intVal = cltiIntegral(Ac,Bc,a,b) OR
% intVal = cltiIntegral(Ac,Bc,a,b,Npoints)
% Calculates the integral int(expm(Act)Bc)dt over the interval [a,b] with
% specified matrices Ac and Bc. Npoints if not specified defaults to the 
% number of steps required so that the step is less than 2*10^-4.
% Otherwise, it is the number of points to be used for the numerical
% integration.
%
% Usage #2:
% intVal = cltiIntegral(contSys,a,b)
% intVal = cltiIntegral(contSys,a,b,Npoints)
% Calculates the integral int(expm(Act)Bc)dt over the interval [a,b] with
% the matrices Ac and Bc included in the system structure contSys so that
% contSys.A=Ac and contSys.B=Bc. Npoints is as before.


%---- Introductory: Check input parameters ----
ni = nargin;
error(nargchk(3,5,ni));
Npoints=-1;
if isstruct(x1)
    error(nargchk(3,4,ni));
    Ac = x1.A;
    Bc = x1.B;
    a=x2;
    b=x3;
    if ni==4
        Npoints=x4;
    end
else
    error(nargchk(4,5,ni));
    Ac=x1;
    Bc=x2;
    a=x3;
    b=x4;
    if ni==5
        Npoints=x5;
    end
end
if Npoints==-1
    Npoints = ceil((b-a)*5000);
end

% Integration step (h):
h=(b-a)/Npoints;


% Integration:
t = a;
intVal=expm(Ac*t)+expm(Ac*b);
for i=1:Npoints-1
    t = t+h;
    intVal = intVal + 2*expm(Ac*t);
end
intVal = intVal*Bc*h/2;