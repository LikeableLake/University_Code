clc; clear; close all;

%% theta method
L=1; %dominio
k=1;
N=30;%numero punti

x=linspace(0,L,N);
h=x(2)-x(1);
hq=h*h;
T=10;   %tempo finale
beta=1/2;
Dt=beta*hq/k;
Nt=round(T/Dt);   %non è preciso perchè ho inserito "troppi dati"

%%
D2= -gallery('tridiag',N-2);
theta=0.5;
I=eye(N-2);
A=I-beta*theta*D2;
B=I+beta*(1-theta)*D2;
Phi0= sin(pi*x/L);
Phi=Phi0(2:N-1)';
for it=1:Nt
    Phi=A\(B*Phi);
    plot(x,[0;Phi;0],'-b',x,Phi0,'r-');drawnow;
end


