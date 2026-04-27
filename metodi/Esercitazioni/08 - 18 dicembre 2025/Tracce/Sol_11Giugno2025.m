close all; clear all; clc;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di risoluzione della traccia del 11 giugno 2025.
% 
% La traccia propone la risoluzione numerica dell'equazione di convezione 
% lineare mediante semidiscretizzazione. La discretizzazione spaziale e' 
% effettuata con uno schema implicito (compatto) e la integrazione 
% temporale con un metodo Runge-Kutta a tre stadi.
%

global Op

L    = 1;                         % Ampiezza del dominio
N    = 50;                        % Numero di punti di discretizzazione
x    = linspace(0,L,N);           % Mesh spaziale
h    = x(2) - x(1);               % Passo spaziale
a    = 1;                         % Velocita' di convezione
C    = 0.8;                       % Numero di Courant
Dt   = C*h/a;                     % Intervallo temporale da Courant
T    = 3*L/a;                     % Tempo finale della simulazione
Nt   = round(T/Dt);               % Numero di step temporali
Dt   = T/Nt;                      % Dt effettivo
Phi0 = exp(sin(2*pi*x/L));        % Condizione iniziale

% Costruzione della matrice di derivazione
%
% Lo schema è dato dalla formula 
%                      A*phi' = B*phi
% per cui la matrice di derivazione è data da
%                      D = inv(A)*B
%
v = zeros(1,N-1);   v(1) = 4;   v(2) = 1;   v(end) = 1;
A = gallery('circul',v); 
v = zeros(1,N-1);               v(2) = 1;   v(end) = -1;
B = gallery('circul',v)*(3/h); 
D = A\B;
% Operatore dello schema
Op = -a*D;
% Ciclo temporale
Phi = Phi0(1:N-1)';
for it = 1:Nt
    t = Dt*it;
% Schema Runge-Kutta
    Phi1 = Phi;                  F1 = RHS(t,Phi1);
    Phi2 = Phi + Dt*F1/3;        F2 = RHS(t,Phi2);
    Phi3 = Phi + Dt*2*F2/3;      F3 = RHS(t,Phi3);
    Phi  = Phi + Dt*(F1 + 3*F3)/4;
% Grafica
    plot(x,Phi0,'r.-'); hold on;
    plot(x,[Phi;Phi(1)],'k.-'); hold off; title('Convezione lineare'); 
    text(.7*L,2,['t = ',num2str(t)]); drawnow;
    IntPhi(it) = sum(Phi)*h;
end
figure(2); 
plot(IntPhi,'.-'); 
title('Integrale della soluzione nel tempo'); drawnow;
