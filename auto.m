function [lambda, vett1, vett2] = auto(A)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[V,D]=eig(A);

lambda=[D(1,1); D(2,2)];
vett1=V(1,:);
vett2=V(2,:);

