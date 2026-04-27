        clear all; close all; clc;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% 15 novembre 2024
% Prof. G. Coppola
% 
% Codice di esempio che risolve l'equazione di Laplace 2D in un dominio
% del tipo 'L' con condizioni al contorno alla Dirichlet.
% La variabile Psi soluzione dell'equazione viene interpretata come la
% funzione di corrente di un moto incompressibile ed irrotazionale, per cui
% le condizioni al contorno sono assegnate in modo da simulare il flusso
% (a potenziale) interno ad un condotto.
% In particolare si prevedono le condizioni al contorno:
%
%                             Psi =1
%                   --------------------------
%    Psi Lineare    |                         |
%     fra 0 e 1     |                         |
%                   |                         |
%                   ---------------           |
%                     Psi = 0      |          | Psi = 1
%                                  |          |
%                                  |          |
%                          Psi = 0 |          |
%                                  |          |
%                                   ----------
%                                   Psi Lineare
%                                    fra 0 e 1
% 
N  = 100;                      % Numero di nodi sul mesh.
L  = 1;                       % Ampiezza del dominio.
x  = linspace(0,L,N); y = x;  % Mesh uniforme su [0,L].
h  = x(2)-x(1);  hq = h*h;    % Passo spaziale (uniforme).

G    = numgrid('L',N);        % Generiamo l'operatore di Laplace...
Lap  = -delsq(G)/hq;          % ...sul dominio a 'L'
Ninc = max(max(G));           % Numero di incognite
p    = zeros(Ninc,1);         % Produzione (nulla)
PsiN = 1;                     % BC sulla parete esterna del condotto 
PSIin  = linspace(1,-1,N);    % BC sulla sezione di inflow
PSIout = linspace(-1,1,N);    % BC sulla sezione di outflow

PSI = nan(size(G));         % Array per la grafica finale
% Assegnazione delle condizioni al contorno
for i = 2:N/2
    for j = 1:N
        if G(i,j)~=0 && G(i-1,j)==0
            p(G(i,j)) = p(G(i,j)) - PsiN/hq;
            PSI(i-1,j) = PsiN;
        end
    end
end
for i = 1:N
    for j = N/2:N-1
        if G(i,j)~=0 && G(i,j+1)==0
            p(G(i,j)) = p(G(i,j)) - PsiN/hq;
            PSI(i,j+1) = PsiN;
        end
    end
end
for i = 1:N
    for j = 1:N/2
        if G(i,j)~=0 && G(i,j-1)==0
            p(G(i,j)) = p(G(i,j)) - PSIin(i)/hq;
            PSI(i,j-1) = PSIin(i);
        end
    end
end
for i = N/2+1:N
    for j = 1:N
        if G(i,j)~=0 && G(i+1,j)==0
            p(G(i,j)) = p(G(i,j)) - PSIout(j)/hq;
            PSI(i+1,j) = PSIout(j);
        end
    end
end
% Risoluzione del sistema
Psi = Lap\p;
% Riversiamo la soluzione nella matrice PSI
PSI(G>0) = Psi;

% Grafica
figure(1); contour(x,y,PSI,20,'w'); axis square
figure(2); surf(x, y, PSI,'EdgeColor','none'); axis square ; colorbar