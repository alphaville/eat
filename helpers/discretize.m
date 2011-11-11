function [Ad Bd] = discretize(Ac,Bc,Ts)
%DISCRETIZE
% Get the discretization of an LTI system
syms x;
Ad = expm(Ac*Ts);
Bd = double(int(expm(Ac*x)*Bc,0,Ts));
clear x;