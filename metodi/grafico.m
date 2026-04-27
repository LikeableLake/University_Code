clc;
clear;
close all;

%%
x=linspace(0, 2*pi, 100);

for i= 1:100
    s=sin(x-i*0.1);
    plot(x,s, '.-w'); pause(0.01)
end;

