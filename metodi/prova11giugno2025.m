clc; close all; clear;
%% definizione del dominio

global f;
L=1;                        %lunghezza dominio
Nx=50;                      %grandezza mesh
x=linspace(0,L,Nx);         %punti x discretizzati
h=x(2)-x(1);  hq=h*h;       %spaziatura mesh spaziale
a=1;                        %velocità convezione
COU=0.8;                    %numero di courant
dt=(COU*h)/a;               %spaziatura mesh temporale
T=3*(L/a);                  %tempo finale
Nt=round(T/dt);             %spaziatura mesh temporale
dt=T/Nt;                    %dt aggiustato per l'arrotondamento
PHI0=exp(sin((2*pi*x)/L));  %condizione iniziale

%% trovo le matrici dello schema
%ipotizzo A*PHI'=(3/h) B*PHI

vA=zeros(1, Nx-1); vA(1)=4; vA(2)=1; vA(end)=1;
vB=zeros(1, Nx-1); vB(2)=1; vB(end)=-1;
A=gallery('circul',vA);
B=gallery('circul',vB);

% PHI'=-a*D*PHI

D=(3/h)*(A\B);

f=-a*D;

PHI=PHI0(1:Nx-1);
PHI=PHI';


%% runge kutta e illustrazione
sumrk4=zeros(Nt,1);
for i=1:Nt
    PHI1=PHI;
    PHI2=PHI+dt*(f*PHI1)/3;
    PHI3=PHI+dt*2*(f*PHI2)/3;
   
    PHI=PHI+dt*(f*PHI1+3*f*PHI3)/4;

    sumrk4(i)=sum(PHI)*h;
    %grafica

    figure(1); subplot(2,1,1) 
    plot(x, PHI0, 'r.-');hold on;
    plot(x, [PHI; PHI(1)], 'w-'); hold off;
    text(0.6,2.2,'convezione lineare con rk4')
    text(0.6, 2, ['ETA: ',num2str(T-(i*dt))]);
    drawnow; 
end


%% metodo ode45
tspan = linspace(0,T,Nt); 
[t, PHI45] = ode45(@odefun, tspan, PHI0(1:Nx-1));

sumODE45=sum(PHI45')'*h;

for i = 1:Nt
    figure (1); subplot(2,1,2); 
    plot(x,PHI0,'.-r'); hold on;
    plot(x, [PHI45(i,:) PHI45(i, 1)], 'w-'); hold off;
    text(0.6,2.2,'convezione lineare con ode45')
    text(0.6, 2, ['ETA: ',num2str(T-(i*dt))]);
    drawnow; 
end

figure (2);subplot (2,1,1);
plot (sumrk4, '.b'); hold on;
subplot (2,1,2);
plot (sumODE45, '.y'); hold on