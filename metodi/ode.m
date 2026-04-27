function [ode] = ode(~,PHI)

A=10*rand(4,4);

ode=A*PHI;
