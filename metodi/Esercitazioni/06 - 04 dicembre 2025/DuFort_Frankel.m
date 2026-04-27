clear; clc; close all

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola 
%
% Codice di esempio che implementa uno schema su tre livelli temporali per
% la risoluzione della equazione di diffusione instazionaria.
% Lo schema (di DuFort e Frankel) ha stencil:
% 
%            t^(n+1)           o     
%                             /|\
%                            / | \
%             t^(n)     o---<  |  >---o
%                            \ | /
%                             \|/
%            t^(n-1)           o
% 
%                    x_(i-1)  x_(i)  x_(i+1)
% 
% e corrisponde al caso n. 8 della tabella di schemi per la diffusione
% tratta dal testo di Richtmyer & Morton (Difference Methods for 
% Initial-Value Problems, 1967, p. 188).
% 
% Lo schema e' dato da:
%
%   (1 + 2*beta)*phi^(n+1) = 2*beta*S*phi^n + (1 - 2*beta)*phi^(n-1),
%
% dove S è la matrice tridiagonale con diagonali [1 0 1].
% Dividendo per (1 + 2*beta) diventa evidente che lo schema è esplicito:
%
%   phi^(n+1) = A*phi^n + B*phi^(n-1)      
%

% Dominio e costanti fisici
L    = 1;                     % Lunghezza del dominio 1D
T    = 1.5;                   % Tempo finale della simulazione
k    = 1;                     % Coefficiente di diffusione
% Discretizzazione
N    = 50;                    % Numero di nodi spaziali
x    = linspace(0,L,N);       % Mesh spaziale
h    = x(2) - x(1);           % Passo spaziale
beta = 20;                    % k*dt/(h^2) 
dt   = beta*(h^2)/k;          % Time step 
Nt   = 150;                   % Numero di istanti temporali 
t    = linspace(0,T,Nt+1);    % Vettore dei tempi
% Eliminiamo primo e ultimo nodo dal mesh (per BCs alla Dirichlet).
x    = x(2:end-1);  x = x(:);
% Condizione iniziale
Phi0 = sin(pi/L*x);
% Operatori differenziali
I = eye(N-2);
S = gallery('tridiag',N-2,1,0,1);
B = 2*beta/(1+2*beta)*S;
C = (1-2*beta)/(1+2*beta)*I;
p = zeros(N-2,1);             % Produzione fisica
q = 2*dt*p;                   % Termine sorgente nello schema
% Condizioni al contorno alla Dirichlet implementate al livello n+1
PhiA   = 0;      PhiB = 0;    % Condizioni al contorno left e right
q(1)   = q(1)  + (2*beta/(1+2*beta))*PhiA;
q(end) = q(end)+ (2*beta/(1+2*beta))*PhiB;
% Costruiamo un array per memorizzare la soluzione nel tempo.
PHI       = zeros(N,Nt+1); % condizione iniziale + Nt time-step 
PHI(:,1)  = [PhiA;Phi0;PhiB];
axis([0 L -0.5 1.5])
% Ciclo instazionario
% Il primo passo temporale lo effettuiamo con un Crank-Nicolson
q_CN     = dt*p;
q_CN(1)  = q_CN(1)   + beta*PhiA;
q_CN(end)= q_CN(end) + beta*PhiB;
D2       = -gallery('tridiag', N-2);
I        = eye(N-2);
A_CN     = I - 0.5*beta*D2;
B_CN     = I + 0.5*beta*D2;
Phi      = A_CN\(B_CN*Phi0 + q_CN);
% Grafica del solo primo step
plot([0;x;L],[PhiA;Phi0;PhiB],'r.-',[0;x;L],[PhiA;Phi;PhiB],'b.-')
title(['it/nt = ',num2str(1),'/',num2str(Nt)])
axis([0 L -.6 1.2]), grid on; drawnow
PHI(:,2) = [PhiA;Phi;PhiB];
% Passi temporali dal secondo in poi
PhiOld = Phi0; 
for it = 2:Nt      
    PhiNew = B*Phi + C*PhiOld + q;     
    PhiOld = Phi;              % il valore 'attuale' diventa 'vecchio',...
    Phi    = PhiNew;           % ...quello 'nuovo' diventa 'attuale'
    PHI(:,it+1) = [PhiA;Phi;PhiB];
% Grafica    
    plot([0;x;L],[PhiA;Phi0;PhiB],'r.-',[0;x;L],[PhiA;Phi;PhiB],'k.-')
    title(['it/nt = ',num2str(it),'/',num2str(Nt)])
    axis([0 L -0.6 1.2]), grid on; drawnow
end
% Grafica tridimensionale
figure; surf(t,[0;x;L],PHI);  shading interp; hold on;
contour3(t,[0;x;L],PHI,50,'k')
xlabel('t'); ylabel('x'); zlabel('\phi')






