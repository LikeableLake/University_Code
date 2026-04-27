clc; close all; clear;

%dati
L=1;
n=60;
x=linspace(0,L,n);
a=1;
k=1;
h=x(2)-x(1);
hq=h^2;

theta=1/2;

D2=-gallery('tridiag', n-2);
D2(end,end)=D2(end, end)+1;
A=gallery("tridiag", n-2, -1, 0, 1);
A(end,end)=-1;
I=eye(n-2, 'like', A);

%% dati che varieranno
i=[0.1 0.3 0.5 0.75 1 1.25 1.5 2 3];
j=1;
for i=i

beta=i;
dt=beta*hq/k;
T=beta*(L^2)/k;
nt=round(T/dt);
dt=T/nt;
t=linspace(0,T,nt);

c=a*dt/h;
omega=2*pi/T;

phiA=cos(omega*t);
neumannB=0;

P=zeros(n-2,nt);
P(1,:)=(k/hq + a/(2*h))*phiA;
P(end,:)=0;

PHI0=cos(2*pi*x/L);
PHI=PHI0(2:n-1)';



for m=1:nt
    PHI=((I/dt - (k*theta/hq)*D2)\(I/dt -a/(2*h)*A + (k*(1-theta)/hq)*D2))*PHI + (I/dt - (k*theta/hq)*D2)\P(:,m);
    figure (1); subplot (5, 5, j); plot(PHI); drawnow; hold off;plot(PHI0); hold on;
end
j=j+1;
end
