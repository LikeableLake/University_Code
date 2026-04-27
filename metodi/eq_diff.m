% discretizzazione delle equazioni differenziali del secondo ordine del
% tipo: -k*d^2phi/dx^2 + a*dphi/dx + b*phi=PI(x)
%la situazione è diversa in base alle condizioni al contorno

close all; clear; clc;
L=1; %definiamo la lunghezza del dominio
a=10;
k=1;
b=1;
PhiA=1; PhiB=1;     %condizioni al contorno dirichlet

%parametri numerici

N=30; %numero punti del dominio
x=linspace(0,L,N); %punti del dominio
h= x(2)-x(1); %il mesh è uniforme quindi h è sempre lo stesso
hq=h*h;

D=gallery("tridiag", N-2, -1, 0, 1)/(2*h); %creo la prima matrice per le derivate prime

D2=gallery("tridiag", N-2, 1,-2,1)/hq; %matrice derivate seconde

%le matrici sono create senza la prima e l'ultima colonna

I=speye(N-2); %matrice identica ma sparsa
p=ones(N-2,1);

A= -k*D2 +a*D +b*I;     %matrice che rappresenta l'equazione differenziale

%se lasciassi così avrei condizione alla dirichlet (produzione pura avrei
%a e b=0 e avrei una parabola)

Phi=A\p;    %vettore soluzione punti interni

plot(x, [0; Phi; 0])        %plotto il grafico di Phi, consierando anche i punti iniziali e finali che in questo caso sono nulli


%% condizioni al contorno

p(1)=p(1) + (k/hq -a/(2*h))*PhiA;
p(N-2)=p(N-2) + (k/hq -a/(2*h))*PhiB;

A= -k*D2 +a*D +b*I;     %matrice che rappresenta l'equazione differenziale

Phi=A\p;
plot(x, [PhiA;Phi;PhiB]);


%attento, se vai troppo preciso il mesh va a puttane

%% condizione alla neumann
g=0; %phi1=phi2-gh

%alla matrice D2 va aggiunto -1 1 alla prima riga
%alla matrice D va aggiunto -1/2h

D(1 ,1)= -1/2*h;
D2(1, 1)=-1/hq;

p(1)=p(1)-(k/h +a/2)*g;

A= -k*D2 +a*D +b*I;     %matrice che rappresenta l'equazione differenziale

Phi=A\p;

plot(x, [Phi(1)-g*h; Phi; PhiB]);


%% mesh variabile


% definisco un h(i) 
%posso definirlo usando una funzuone e decido io dove è più denso

csi = linspace(0,L,N);
x = L* (csi/L).^2;  %infittisce vicino a 0
x = L* (csi/L).^0.5; %infittisce vicino a L

%chiamo la pesider

for i=3:N-2
    ir=i-1;
    xs = x(i-1:i+1);
    xc =xs(2);
    w1= PesiDer(xs,xc,1);   D(ir, ir-1:ir+1)=w1;
    w2= PesiDer(xs,xc,2);   D2(ir, ir-1:ir+1)=w2;
end

%il ciclo fot riempie l'interno della matrice, la prima e ultima componente
%le faccio a parte


