clc; close all; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
%  
% Codice di esempio che simula l'equazione di diffusione instazionaria in 
% 2 dimensioni (2D) in un quadrato con varie condizioni al contorno alla 
% Dirichelet o Neumann-Dirichlet. 
% Per la integrazione temporale si utilizza uno schema 'theta-method'.
% La equazione è:
% 
%                   phi_t = k*(phi_xx + phi_yy) + p
% 
% La codifica generale dello schema e' data da:
% 
%                  A*phi(n+1) = B*phi(n) + q
% 
% Il codice implementa due possibili boundary conditions:
% 1-2) Dirichlet BCs (anche non omogenee) su tutti i lati del dominio
% 3) Dirichlet su due lati del dominio e Neumann sugli altri due.
% 
L    = 1;                     % Estensione del dominio (lato del quadrato)
N    = 40;                    % Numero di punti in ciascuna direzione
k    = 1;                     % Diffusivita'
x    = linspace(0,L,N);       % Mesh nella direzione x
y    = x;                     % Mesh nella direzione y
h    = x(2)-x(1);             % Passo spaziale
G    = numgrid('S',N);        % Generazione della griglia nel quadrato
Lap  = -delsq(G);             % Operatore di Laplace con BCs alla Dirichlet
beta = 0.8; teta = 0.5;       % Parametri numerici
Dt   = beta*h^2/k;            % Passo temporale
P    = 5*ones((N-2)^2,1);     % Termine sorgente nella equazione
P    = P*Dt;                  % Termine sorgente nella discretizzazione
Phi0 = ones((N-2)*(N-2),1);   % Condizione iniziale
Case = 1;                     % Selezione del caso
switch Case
    case 1
        % Condizione al contorno (Dirichlet)
        D_West = 0;            D_Est = 0;
        D_Nord = 0;            D_Sud = 0;
        % Condizione al contorno (Dirichlet)
        for k = 1:N-2
            iWest =  k;                P(iWest) = P(iWest) + beta*D_West;
            iEst  = (N-2)*(N-3)+k;     P(iEst)  = P(iEst)  + beta*D_Est;
            iNord = (k-1)*(N-2)+1;     P(iNord) = P(iNord) + beta*D_Nord;
            iSud  =  k*(N-2);          P(iSud)  = P(iSud)  + beta*D_Sud;
        end
    case 2
        % Condizione al contorno (Dirichlet)
        D_West = 0;            D_Est = 0;
        D_Nord = 1;            D_Sud = 1;
        % Condizione al contorno (Dirichlet)
        for k = 1:N-2
            iWest =  k;                P(iWest) = P(iWest) + beta*D_West;
            iEst  = (N-2)*(N-3)+k;     P(iEst)  = P(iEst)  + beta*D_Est;
            iNord = (k-1)*(N-2)+1;     P(iNord) = P(iNord) + beta*D_Nord;
            iSud  =  k*(N-2);          P(iSud)  = P(iSud)  + beta*D_Sud;
        end
    case 3
        % Condizione al contorno (N: Neumann, D: Dirichlet)        
        N_Nord = -1;           N_Sud = 0;
        D_West = 1;            D_Est = 0;
        for k = 1:N-2
        % Correzione termine sorgente per le BC    
            iWest =  k;                P(iWest) = P(iWest) + beta*D_West;            
            iEst  = (N-2)^2 - k + 1;   P(iEst)  = P(iEst)  + beta*D_Est;                        
            iNord = (k-1)*(N-2)+1;     P(iNord) = P(iNord) + N_Nord*h*beta;
            iSud  =  k*(N-2);          P(iSud)  = P(iSud)  + N_Sud*h*beta;
        % Correzione matrice per le BC (solo nel caso Neumann)
            Lap(iNord,iNord) = Lap(iNord,iNord) + 1;
            Lap(iSud,iSud)   = Lap(iSud,iSud)   + 1;
        end
end
% Definizione degli operatori per la discretizzazione spaziale
I    = eye(size(Lap));          % Matrice identica
A    = I - beta*teta*Lap;       % Operatore del 'theta-method' al LHS
B    = I + beta*(1-teta)*Lap;   % Operatore del 'theta-method' al RHS
% Ciclo instazionario
Phi  = Phi0;
for it=1:800
    Phi = A\(B*Phi + P);            % Soluzione del sistema implicito
    PHI = reshape(Phi,N-2,N-2);     % Reshape per la grafica
    switch Case                     % Aggiunta delle BCs
        case {1,2}
            bcWest  = D_West*ones(N-2,1);      bcEst = D_Est*ones(N-2,1);
            bcNord  = D_Nord*ones(1,N);        bcSud = D_Sud*ones(1,N);
            PHI     = [bcWest,PHI,bcEst];
            PHI     = [bcNord;PHI;bcSud];            
        case 3
            bcNorth = PHI(1,:)   + N_Nord*h;
            bcSouth = PHI(end,:) + N_Sud*h;            
            PHI     = [bcNorth;PHI;bcSouth];
            bcWest  = D_West*ones(N,1);
            bcEst   = D_Est*ones(N,1);            
            PHI     = [bcWest,PHI,bcEst];
    end
    % Grafica
    figure(1); surfl(x,y,PHI); shading interp;
    hold on; contour3(x,y,PHI,'k');  
    axis([0 L 0 L 0 3])
    drawnow;  hold off
end