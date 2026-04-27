function [odefunPHIE] = odefunPHIE(~,PHIE)

f=[-2 -1; 1 0];

odefunPHIE=[1;0]-f*PHIE;