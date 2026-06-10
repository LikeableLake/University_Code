function [eq] = par(A,B,C,D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    eq=1/((1/A)+(1/B));
elseif nargin==3
    eq=1/((1/A)+(1/B)+(1/C));
elseif nargin==4
    eq=1/((1/A)+(1/B)+(1/C)+(1/D));
end



