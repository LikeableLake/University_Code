clc; close all; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che simula la equazione di convezione lineare con
% condizioni al contorno periodiche.
% Il codice utilizza un approccio semidiscretizzato, nel quale l'equazione
% viene prima discretizzata nello spazio e poi integrata nel tempo. In
% questo script si implementa una classica procedura RK4
%  

global A               % Dichiarazioni globali

L  = 1;                 % Lunghezza del dominio
N  = 60;                % Numero di punti di discratizzazione spaziale
x  = linspace(0,L,N);   % Definizione del mesh
h  = x(2) - x(1);       % Passo spaziale
a  = 1;                 % Velocità di convezione
C  = 1.5;               % Numero di Courant
T  = L/a;               % Tempo finale
Dt = h*C/a;             % Dt per la stabilita'
Nt = round(T/Dt);       % Numero di passi temporali
Dt = T/Nt;  C = a*Dt/h; % Rideterminazione del Dt e C

% Matrice di derivata prima su mesh periodico
v  = zeros(1,N-1);   v(2) = 1;   v(end) = -1;
D  = gallery('circul',v)/(2*h);
A = -a*D;
% Condizione iniziale
x0  = L/2;              Phi0 = exp(-100*(x-x0).^2);   
Phi = Phi0(1:end-1)'; 
% Procedura RK4
for it = 1:Nt
    t = it*Dt;
    Phi1 = Phi;                   F1 = ConvFun(t,Phi1);
    Phi2 = Phi + Dt*0.5*F1;       F2 = ConvFun(t,Phi2);
    Phi3 = Phi + Dt*0.5*F2;       F3 = ConvFun(t,Phi3);
    Phi4 = Phi + Dt*F3;           F4 = ConvFun(t,Phi4);

    Phi = Phi + Dt*( (1/6)*(F1+F4) + (1/3)*(F2+F3) );
% Grafica
    plot(x,Phi0,'.-r',x,[Phi; Phi(1)],'.-k'); 
    axis([0 L -0.5 1.5]); drawnow;
end
