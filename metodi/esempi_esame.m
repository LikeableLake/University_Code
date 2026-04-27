clc; close all; clear;

%%

%comandi di matlab che discretizzano l'operatore laplaciano per domini noti

N=32;  %numero punti
L=1; %lato del quadrato

x=linspace (0, L, N);
h=x(2)-x(1);
G= numgrid('S', N); %crea una griglia con condizioni dirichlet omogenee quadrata

L= -delsq(G);   %operatore del laplaciano per quadrato

%ora definisco le condizioni al contorno

Ninc = max(max(G)); %il numero di incognite è il numero massimo nella matrice, si fa due volte perchè la primna volta dà il vettore massimo
q=zeros(N);
PhiB = 1; %valore dirichlet di bordo

%dobbiamo modificare l'operatore di laplace per avere le condizioni al
%contorno giuste

for i= 1:N
    for j= 1:N
        if(G(i,j)~=0 && G(i,j+1)==0)
            q(i,j)= q(i,j)-PhiB/(h^2);
        end
    end
end

q= q(2:end-1, 2:end-1);
q=q(:);

%ordino q lessicografico


Phi= L\q;
Phi = reshape(Phi, N-2, N-2);   %riformo la funzione a quadrato solo punti interni

%per la funzine completa necessito anche dei valori di contorno che ho
%detto essere 1 a destra e 0 negli altri tre lati

Phi= [zeros(N-2, 1) Phi ones(N-2, 1)*PhiB];

Phi=[zeros(1, N); Phi; zeros(1,N)];

%surf(Phi); per visualizzare la funzione

%condizioni al contorno diverse 

%% seconda parte

q=Phi(2:(end-1), 2:(end-1));
q=q(:);
Psi = L\q;
Psi = reshape(Psi, N-2, N-2);

%% terza parte dominio rettangolare

%spy() ti fa vedere la matrice in un grafico più caruccia 


H=numgrid('H', 8);
spy(H)

%sicurop esce all'esame condizioni non alla dirichlet
