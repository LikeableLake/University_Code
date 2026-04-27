close all; clear all; clc;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice che risolve la equazione di convezione-diffusione-reazione 
% stazionaria monodimensionale con schema assegnato in modo arbitrario.
% Si prevede la discretizzazione su un mesh eventualmente non uniforme. 
% Nella codifica originale si prevedono i seguenti parametri di
% discretizzazione:
% 1) il mesh x e' non uniforme e lo si ottiene equispaziando la variabile
%    csi = x^p (x = csi^(1/p))
% 2) Per la diffusione e la convezione si utilizza uno schema centrale su 
%    tre punti in modo che le derivate nel nodo i-mo siano calcolate
%    utilizzando i valori sui nodi i-1, i e i+1
% 5) Le condizioni al contorno sono assegnate alla Dirichlet.
% 
% Il codice è facilmente modificabile per prevedere diversi schemi di
% discretizzazione. Si suggerisce ad esempio di implementare: schemi upwind
% per la convezione, condizioni al contorno diverse (e.g. Neumann o Robin),
% mesh variabili con leggi diverse, etc.
% 
% Nella scrittura del codice si fa riferimento allo schema per il mesh:
% 
%                        x = 0                           x = L
%                          o-------x------x-- ... --x---x--o
% Numerazione nodi:        1       2      3        N-2 N-1 N
% 
% Il codice utilizza la function PesiDer.m per il calcolo dei pesi delle
% formule di derivazione.
% 

% Parametri fisici
L    = 1;                        % Lunghezza del dominio.
k    = 0.1;                      % Coefficiente di diffusione.
a    = 5;                        % Coefficiente di convezione.
b    = 1;                        % Coefficiente di reazione.
PhiA = 1;   PhiB  = 2;           % Boundary conditions.
% Parametri numerici
N    = 30;                       % Numero di nodi di discretizzazione.
csi  = linspace(0,L,N);          % Mesh csi (uniforme) 
x    = L*(csi/L).^0.5;           % Mesh x di discretizzazione 
% Grafica del mesh.
figure;   plot(x,0*x,'k.-','markersize',12);   title('Mesh')
%% Inizializzazione delle matrici di derivazione prima e seconda.
D1   = sparse(N-2);      D2 = sparse(N-2);           
I    = speye(N-2);               % Matrice identica.
p    = zeros(N-2,1);             % Produzione.
% Calcolo delle componenti delle matrici 
% Primo nodo interno.
in   = 2;                        % Indice di nodo 
ir   = in - 1;                   % Indice di riga nelle matrici.
xs   = x(in-1:in+1);             % Stencil per la derivata prima.
xc   = xs(2);                    % Collocazione per la derivata prima.
w1   = PesiDer(xs,xc,1);         % Pesi della derivata prima.
D1(ir,ir:ir+1) = w1(2:3);        % Riempimento della matrice.
w10  = w1(1);                    % Conserviamo il primo 'peso' per la BC
xs   = x(in-1:in+1);             % Stencil per la derivata seconda.
xc   = xs(2);                    % Collocazione per la derivata seconda.
w2   = PesiDer(xs,xc,2);         % Pesi della derivata seconda.
D2(ir,ir:ir+1) = w2(2:3);        % Riempimento della matrice.
w20  = w2(1);                    % Conserviamo il primo 'peso' per la BC
% Calcolo delle componenti delle matrici per i punti interni.
for in = 3:N-2
    ir = in - 1;                 % Indice di riga nelle matrici.
    xs = x(in-1:in+1);           % Stencil per la derivata prima.
    xc = xs(2);                  % Collocazione per la derivata prima.
    w1 = PesiDer(xs,xc,1);       % Pesi della derivata prima.
    D1(ir,ir-1:ir+1) = w1(1:3);  % Riempimento della matrice.
    xs  = x(in-1:in+1);          % Stencil per la derivata seconda.
    xc  = xs(2);                 % Collocazione per la derivata seconda.
    w2  = PesiDer(xs,xc,2);      % Pesi della derivata seconda.
    D2(ir,ir-1:ir+1) = w2(1:3);  % Riempimento della matrice.
end
% Ultimo nodo interno.
in   = N-1;                      % Indice di nodo 
ir   = in - 1;                   % Indice di riga nelle matrici.
xs   = x(in-1:in+1);             % Stencil per la derivata prima.
xc   = xs(2);                    % Collocazione per la derivata prima.
w1   = PesiDer(xs,xc,1);         % Pesi della derivata prima.
D1(ir,ir-1:ir) = w1(1:2);        % Riempimento della matrice.
w1L  = w1(3);                    % Conserviamo il primo 'peso' per la BC
xs   = x(in-1:in+1);             % Stencil per la derivata seconda.
xc   = xs(2);                    % Collocazione per la derivata seconda.
w2   = PesiDer(xs,xc,2);         % Pesi della derivata seconda.
D2(ir,ir-1:ir) = w2(1:2);        % Riempimento della matrice.
w2L  = w2(3);                    % Conserviamo il primo 'peso' per la BC
% Operatore convezione-diffusione-reazione
A    = -k*D2 + a*D1 + b*I;     
% Boundary conditions
p(1)   = p(1)   + (k*w20 - a*w10)*PhiA;
p(N-2) = p(N-2) + (k*w2L - a*w1L)*PhiB;
% Risoluzione del sistema.
Phi  = A\p;                       
% Grafica della soluzione.
figure; plot(x,[PhiA;Phi;PhiB],'w.-'); axis([0 L .5 2.5])       
xlabel('x'); ylabel('\phi');           title('Grafica della soluzione')
Pe = a*(max(diff(x)))/(2*k);     % Calcolo del massimo Peclet di cella.
text(0.1,1.2,['Peclet di cella = ',num2str(Pe)]);



