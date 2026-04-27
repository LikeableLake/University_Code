clear all; clc; close all      

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice che implementa la soluzione di un problema di diffusione-convezione
% stazionaria 1D con produzione assegnata p(x):
% 
%                        -k*phi'' + a*phi' = p(x)
%                      
% mediante schemi numerici del secondo ordine su mesh uniforme con 
% condizioni al contorno alla Dirichlet-Dirichlet e Neumann-Dirichlet.
% 
% 
L = 1;                            % Ampiezza del dominio.
N = 30;                           % Numero di nodi sul mesh.
x = linspace(0,L,N);              % Mesh uniforme su [0,L].
h = x(2) - x(1);      hq = h*h;   % Passo spaziale (uniforme).
k = 1;                            % Coefficiente di diffusività.
a = 10;                           % Velocita' di convezione.
% 
% Per il mesh si utilizza la notazione:
% 
%                                                   |- h -| 
%                           x = 0                   |     |   x = L
%                             o-----x-----x-- ... --x-----x-----o
% Numerazione nodi:           1     2     3        N-2    N-1   N
% Numerazione incognite:            1     2        N-3    N-2   
% 
% I nodi contrassegnati dal simbolo 'o' sono nodi su cui si assegna la
% condizione al contorno.
% 
%--------------------------------------------------------------------------
% Caso con BCs Dirichlet-Dirichlet
PhiA = 1;                         % BC a sinistra. 
PhiB = 2;                         % BC a destra. 
% Matrici di derivazione
D1 = gallery('tridiag',N-2,-1,0,1)/(2*h);
D2 = gallery('tridiag',N-2,1,-2,1)/hq;
% Costruzione dell'operatore.
A  = -k*D2 + a*D1;
p  = ones(N-2,1);                 % Assegno la produzione
% Il termine noto deve tenere in conto anche delle condizioni al contorno
% (non omogenee in generale) date alla Dirichlet.
p(1)   = p(1)   + (k/hq + a/(2*h))*PhiA;
p(end) = p(end) + (k/hq - a/(2*h))*PhiB;
% Risoluzione
Phi  = A\p;
Phi  = [PhiA; Phi; PhiB];         % 'Orlo' la soluzione con le BCs.
% Grafica 
subplot(2,1,1); plot(x,Phi,'k.-');  
xlabel('x');       ylabel('\phi(x)')
title(['BCs Dirichlet-Dirichlet']);
%--------------------------------------------------------------------------
% Caso con BCs Dirichlet-Neumann
PhiA = 1;                         % BC a sinistra (Dirichlet). 
g    = 1;                         % BC a destra (Neumann). 
%
% NB la BC alla Neumann in x=0 e' data nella forma: 
%                    -k*phi' = g
% e discretizzata con una formula forward
%
% Matrici di derivazione
D1 = gallery('tridiag',N-2,-1,0,1)/(2*h);
D2 = gallery('tridiag',N-2,1,-2,1)/hq;
% Modifica delle matrici di derivazione per effetto delle BCs
D1(1,1) = -1/(2*h);       D2(1,1) = -1/hq; 
% Costruzione dell'operatore.
A  = -k*D2 + a*D1;
p  = zeros(N-2,1);                 % Assegno la produzione
% Il termine noto deve tenere in conto anche delle condizioni al contorno
% (non omogenee in generale) 
p(1)   = p(1)   + (1/h  + a/(2*k))*g;
p(end) = p(end) + (k/hq - a/(2*h))*PhiB;
% Risoluzione
Phi  = A\p;
Phi  = [Phi(1) + g*h/k; Phi; PhiB]; % 'Orlo' la soluzione con le BCs.
% Grafica 
subplot(2,1,2); plot(x,Phi,'k.-');  
xlabel('x');       ylabel('\phi(x)')
title(['BCs Neumann-Dirichlet']);

