clc; clear; close all;

k=1;

% discretizzazione dominio spaziale

L=1;
N=60;
x=linspace(0,L,N);
h=x(2)-x(1);
hq=h^2;

% discretizzazione dominio temporale

beta=0.4;
dt=beta*hq/k;
T=1;
Nt=round(T/dt);
t=linspace(0,T,Nt);
dt=t(2)-t(1);

% produzione

f=sin(pi*x/L)';
%f è costante nel tempo

% condizione iniziale

PHI0=cos(2*pi*x/L);
PHI=PHI0(2:N-1)';

% matrice dell'equazione e condizioni al contorno

PHIA=exp(-k*t)';         %dirichlet a x=0
gNEUMANN=0;             %neumann nulla a x=L

D2=-gallery('tridiag', N-2);

%per dirichlet non necessito di cambiare le condizioni al contorno, per
%neumann devo usare lo schema a tre punti quindi:

D2(end, end-1)=2/3;     D2(end,end)=-2/3;

%matrici del theta method

theta=0.5;
I=eye(N-2);
A=I+beta*(1-theta)*D2;
B=I-beta*(theta)*D2;
Dir=zeros(N-2,1); Dir(1)=1;
Neu=zeros(N-2,1); Neu(end)=1;

%soluzione

PHI=B\(A*PHI+dt*f(2:N-1)+Neu*gNEUMANN*dt/h+Dir*beta*(theta*PHIA(1)+(1-theta)*PHIA(2)));
figure(2); plot(x(2:N-1), PHI); drawnow;


% %%
% for it=2:Nt
%     PHI=B\(A*PHI+dt*f(2:N-1)+Neu*gNEUMANN*dt/h+Dir*beta*(theta*PHIA(it-1)+(1-theta)*PHIA(it)));
%     figure(2); plot(x, [PHIA(it); PHI; PHI(end)]); drawnow; hold off;
%     plot(x, PHI0); drawnow; hold on;
% end

%% condizione alla robin

hROBIN=exp(-k*t)';

% i coefficienti sulla prima riga di D2 saranno:

D2(1, 1)=(-4*h+2)/(2*h-3);     D2(1,2)=(-2+2*h)/(2*h-3);

A=I+beta*(1-theta)*D2;
B=I-beta*(theta)*D2;
Rob=zeros(N-2,1); Rob(1)=1;
Neu=zeros(N-2,1); Neu(end)=1;

%soluzione

PHIr=PHI0(2:N-1)';
PHIr=B\(A*PHIr+dt*f(2:N-1)+Neu*gNEUMANN*dt/h+Rob*beta*(theta*hROBIN(1)+(1-theta)*hROBIN(2)));
figure(3); plot(x(2:N-1), PHIr); drawnow;

%%

for it=2:Nt
    PHIr=B\(A*PHIr+dt*f(2:N-1)+Neu*gNEUMANN*dt/h+Rob*beta*(theta*hROBIN(it-1)+(1-theta)*hROBIN(it))/(1-3/(2*h)));
    figure(3); plot(x, [PHIr(1); PHIr; PHIr(end)]); drawnow; hold off;
    plot(x, PHI0); drawnow; hold on;
end