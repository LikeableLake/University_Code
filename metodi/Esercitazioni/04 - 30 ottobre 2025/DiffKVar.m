clc; close all; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice che risolve la equazione di diffusione stazionaria 
% monodimensionale con coefficiente di diffusione variabile.
% 
% La equazione e' posta nella forma:
% 
%      - d(k(x)dPhi/dx)/dx = p 
% 
% e le condizioni al contorno sono assegnate alla Dirichlet.
% 

global L                      % Dichiarazione variabili globali.

N    = 40;                    % Numero di nodi sul mesh.
L    = 1;                     % Ampiezza del dominio.
x    = linspace(0,L,N);       % Mesh uniforme su [0,L].
h    = x(2)-x(1);  hq = h*h;  % Passo spaziale (uniforme).
p    = zeros(N-2,1);          % Produzione.  
PhiA = 1;                     % BC in x=0
PhiB = 2;                     % BC in x=L 

% Assegnazione sparsa della matrice della discretizzazione.
DkD  = sparse(N-2);
% Per la discretizzazione dell'equazione si utilizza lo schema:
% 
%                               i-1/2   i+1/2      |-- h --| 
%    x = 0                       |       |         |       |     x = L
%      o---|---x---|---x ... x---|---x---|---x ... x---|---x---|---o
%      1       2       3    i-1      i      i+1   N-2     N-1      N
% 
% dove i valori ku' vengono stimati nei nodi i-1/2 e i+1/2 con formule di
% derivazione su due punti di ampiezza h. Successivamente, i valori della
% derivata (ku')' vengono stimati di nuovo nei nodi i-mi mediante utilizzo
% ancora di formule di derivazione su due punti di ampiezza h basate sui
% precedenti valori di ku' nei nodi i-1/2 e i+1/2.
% 
% Per la stima dei valori di ku' nei nodi i-1/2 e i+1/2 c'e' bisogno di una
% function 'kappa' che calcola i valori della funzione k(x) in una
% arbitraria posizione del mesh, anche non nodale, a partire da una legge
% analitica.

in = 2;                        % Primo nodo interno
ir = in-1;                     % Indice della riga della matrice
km = kappa(x(in)-h/2);         % Valore di k in i-1/2.
kp = kappa(x(in)+h/2);         % Valore di k in i+1/2.
% Pesi per lo schema numerico.
w0         = -km/hq;           % Da conservare per le BC
DkD(ir,ir) = (km+kp)/hq;    DkD(ir,ir+1) = -kp/hq;   
% Riempimento dei pesi della matrice nei punti interni. 
for in = 3:N-2
    ir = in-1;                 % Indice della riga della matrice
    km = kappa(x(in)-h/2);     % Valore di k in i-1/2.
    kp = kappa(x(in)+h/2);     % Valore di k in i+1/2.
% Pesi per lo schema numerico.
    DkD(ir,ir-1) = -km/hq;  DkD(ir,ir) = (km+kp)/hq;   DkD(ir,ir+1) = -kp/hq;   
end
in  = N-1;                     % Ultimo nodo interno
ir = in-1;                     % Indice della riga della matrice
km = kappa(x(in)-h/2);         % Valore di k in i-1/2.
kp = kappa(x(in)+h/2);         % Valore di k in i+1/2.
% Pesi per lo schema numerico.
DkD(ir,ir-1) = -km/hq;      DkD(ir,ir) = (km+kp)/hq;
wL           = -kp/hq;         % Da conservare per le BC

% Assegnazione delle condizioni al contorno.
p(1)       = p(1)   - w0*PhiA;
p(end)     = p(end) - wL*PhiB;

% Risoluzione del sistema
Phi=DkD\p;

% Grafica.
figure; plot(x,[PhiA;Phi;PhiB],'k.-')  
xlabel('x'); ylabel('\phi'); title('Grafica della soluzione')

