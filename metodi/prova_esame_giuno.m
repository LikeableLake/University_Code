close all; clc; clear;

L=1;
N=50;
x = linspace(0, L, N); h=x(2)-x(1);
c=0.8;
a=1;
Dt=c*h/a;
T=L/a;
Nt=round(T/Dt); Dt=T/Nt; c=a*Dt/h;

Phi0= exp(sin(2*pi*x/L));


v= zeros(1,N-1); v(1)=4; v(2)= 1; v(end)=1;
A=gallery('circul', v);

u= zeros(1,N-1); u(2)=1; u(end)= -1;
B=gallery('circul', v);

D=A\B;
Phi=Phi0(1:N-1)';
for it =1:Nt
    Phi1=Phi;
    Phi2=Phi+Dt*F1/3;
    Phi3=Phi+Dt*F2*2/3

    Phi= Phi+ Dt*(F1 +3*F3)

    %la metterà su teams...