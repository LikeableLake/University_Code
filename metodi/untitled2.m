close all; clear; clc
%utilizzo di pesider
h=0.1;
x=[-h; 0; h];
f=exp(x);

w= PesiDer(x,x(2),2);

dfn=w*f; %derivata seconda numerica

%doc logspace
%esercizio sul trovare l'errore tra derivata numerica e analitica
