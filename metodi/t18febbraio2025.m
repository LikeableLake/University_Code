clc; close all; clear;

%% dati rettangolo
AR=3;
Lx=3;
Ly=Lx/AR;

Nx=(AR*30)+1;                       %modificare il numero moltiplicato per AR 
Ny=1+(Nx-1)/AR;                     %per ottenere Ny e Nx interi con hx=hy

x=linspace(0,Lx,Nx);
y=linspace(0,Ly,Ny);

h=x(2)-x(1);
hq=h^2;

%% costruzione rettangolo
% il rettangolo deve avere Lx righe e Ly colonne quindi

G=1:(Nx-2)*(Ny-2);
G=reshape(G, Ny-2, Nx-2);

G=[zeros(Ny-2,1) G zeros(Ny-2,1)];
G=[zeros(1,Nx); G ; zeros(1,Nx)];

%% costruzione laplaciano

L=-delsq(G);

%% costruzione funzione di sorgente

P=zeros(Ny-2,Nx-2);

for i=1:Nx-2
    for j=1:Ny-2
        P(j,i)=-(y(j+1)-Ly);
    end
end

%% condizioni al contorno per il laplaciano

for i=1:Ny
    for j=1:Nx
        node=G(i,j);
        if (node~=0 && G(i-1, j)==0)
            L(node,node)=L(node,node)+1;
        end
    end
end

%% condizioni al contorno per la sorgente

P(:,1)= P(:,1) + 0;         % dalle condizioni alla neumann a nord
P(:,end)= P(:, end) + 0;   % dalle dirichlet a sud

P(1,:)= P(1, :) - (y(1)*(y(1)-2*Ly))/hq;               %est
P(end, :)= P(end, :) - (y(end)*(y(end)-2*Ly))/hq;       %ovest

%% risoluzione

p=P(:);
L=L/hq;

phi=L\p;

PHI=reshape(phi, Ny-2, Nx-2);



figure(1); surf( PHI', 'FaceColor', 'interp');
colorbar; title('soluzione all''equazione di poisson'); xlabel('x');
ylabel('y'); zlabel('\Phi', 'Rotation', 0);

figure(2); contour( PHI'); title('contour lines'); xlabel('x');
ylabel('y', 'Rotation', 0);

%% seconda parte

%soluzione dell'equazione differenziale
%l'equazione diventa un sistema di due equazioni:

%phie'' +2phie' +phie = 1
%psie = phie', psie' = phie''

%psie'+2psie +phie = 1
%psie = phie'

%psie' = 1 -2psie -phie
%phie' = psie

% [psie';phie']= [1;0] - [-2 -1 ; 1 0]*[psie;phie]

% chiamerò f la funzione matriciale e PHIE il vettore [psie;phie]

% devo integrare su un timespan 'spaziale' che è y in quanto PHIE è
% funzione di solo y

%% rk4

f=[-2 -1; 1 0];
I=[1; 0];

PHIE=NaN(2, Ny);
PHIE0=[1;0];
PHIE(:,1)=PHIE0;

for n=1:Ny-1

PHIE1=PHIE(:,n);
PHIE2=PHIE(:,n)+h*(I-f*PHIE1)/2;
PHIE3=PHIE(:,n)+h*(I-f*PHIE2)/2;
PHIE4=PHIE(:,n)+h*(I-f*PHIE3);

PHIE(:,n+1)=PHIE(:,n)+h*(I-(f*PHIE1/6 +f*PHIE2/3 +f*PHIE3/3+f*PHIE4/6));
end

PHIE=PHIE(2,:)';

P(:, 1)= P(:, 1) - PHIE(2:end-1)/hq;   %nuova condizione al contorno

p=P(:);
L=L/hq;

phi=L\p;

PHI=reshape(phi, Ny-2, Nx-2);



figure(3); surf( PHI', 'FaceColor', 'interp');
colorbar; title('soluzione all''equazione di poisson'); xlabel('x');
ylabel('y'); zlabel('\Phi', 'Rotation', 0);

figure(4); contour( PHI'); title('contour lines'); xlabel('x');
ylabel('y', 'Rotation', 0);
