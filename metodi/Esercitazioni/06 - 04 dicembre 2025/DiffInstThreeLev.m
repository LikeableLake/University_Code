clc, clear all, close all

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che implementa uno schema su tre livelli temporali per
% la risoluzione della equazione di diffusione instazionaria.
% Lo schema ha stencil:
% 
%            t^(n+1)   o----o----o
%                           |
%                           |
%             t^(n)         o     1 + theta
%                           |
%                           |
%            t^(n-1)        o      - theta
% 
% 
% e corrisponde al caso n. 10 della tabella di schemi per la diffusione
% distribuita a lezione, tratta dal testo di Richtmyer & Morton (Difference 
% Methods for Initial-Value Problems, 1967, p. 188).
% 
% Lo schema ha codifica:
% 
%                  A*f(n+1) = B*f(n) + C*f(n-1) + q
% 
% e risulta essere incondizionatamente stabile.
% 

% Dominio e costanti fisici
L    = 1;                 % Lunghezza del dominio 1D
T    = .5;                % Tempo finale della simulazione
k    = 1;                 % Coefficiente di diffusione
% Parametri della iscretizzazione
N    = 50;                % Numero di nodi spaziali
beta = 5;                 % k*dt/(dx^2) e da questo si ricava il dt
x    = linspace(0,L,N);   % Mesh spaziale
h    = x(2) - x(1);       % Passo spaziale
dt   = beta*(h^2)/k;      % time step 
Nt   = round(T/dt);       % Arrotondo T/dt, il tempo finale potrebbe... 
                          % ... essere leggermente diverso da T
theta = 0.25;             % Parametro dello schema.
% Eliminiamo primo e ultimo nodo dal mesh (per BCs alla Dirichlet).
x    = x(2:end-1);  x = x(:);
% Condizione iniziale
Phi0 = cos((3/2)*pi*x/L);% + .5*sin(3*pi*x/L);
% Operatori differenziali
I    = eye(N-2);                    % Matrice identità
D2   = -gallery('tridiag',N-2);     % Matrice di derivazione seconda
A    = (1+theta)*I-beta*D2;         % Matrice A dello schema
B    = (1+2*theta)*I;               % Matrice B ...
C    = -theta*I;                    % Matrice C ...
p    = zeros(N-2,1);                % Produzione fisica
q    = dt*p;                        % Termine sorgente nello schema
% Condizioni al contorno alla Dirichlet implementate al livello n+1
PhiA   = 1;      PhiB = 0;          % Condizioni al contorno left e right
q(1)   = q(1)   + beta*PhiA;        % Modifica del termine noto per...
q(end) = q(end) + beta*PhiB;        % ... tenere conto delle BCs
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
plot([0;x;L],[PhiA;Phi0;PhiB],'r.-',[0;x;L],[PhiA;Phi;PhiB],'k.-')
title(['it/nt = ',num2str(1),'/',num2str(Nt)])
axis([0 L -1.2 1.2]), grid on; drawnow
% Passi temporali dal secondo in poi
PhiOld = Phi0;
for it = 2:Nt
    PhiNew = A\(B*Phi + C*PhiOld + q);
    PhiOld = Phi;              % il valore 'attuale' diventa 'vecchio',...
    Phi    = PhiNew;           % ...quello 'nuovo' diventa 'attuale'
% Grafica    
    plot([0;x;L],[PhiA;Phi0;PhiB],'r.-',[0;x;L],[PhiA;Phi;PhiB],'k.-')
    title(['it/nt = ',num2str(it),'/',num2str(Nt)])
    axis([0 L -1.2 1.2]), grid on; drawnow
end
