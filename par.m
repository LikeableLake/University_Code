function [eq] = par(A,B,C)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin==2
    eq=1/((1/A)+(1/B));
else if nargin==3
    eq=1/((1/A)+(1/B)+(1/C));
end



