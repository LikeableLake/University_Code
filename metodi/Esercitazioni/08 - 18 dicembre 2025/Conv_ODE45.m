clc; close all; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che simula la equazione di convezione lineare con
% condizioni al contorno periodiche.
% Il codice utilizza un approccio semidiscretizzato, nel quale l'equazione
% viene prima discretizzata nello spazio e poi integrata nel tempo. In
% questo script si utilizza il solver ODE45 predefinito in MATLAB
% 

global A                % Dichiarazioni globali

L  = 1;                 % Lunghezza del dominio
N  = 60;                % Numero di punti di discratizzazione spaziale
x  = linspace(0,L,N);   % Definizione del mesh
x  = x(1:end-1);        % Condizioni al contorno periodiche...
h  = x(2) - x(1);       % Passo spaziale
a  = 1;                 % Velocità di convezione
T  = L/a;               % Tempo finale
Nt = 100;               % Numero di passi temporali per l'output di ODE45
% Matrice di derivata prima su mesh periodico
v = zeros(1,N-1); v(2) = 1; v(end) = -1;
D = gallery('circul',v)/(2*h);
A = -a*D;
% Parametri per il solver ODE45
Phi0 = exp(-100*(x-L/2).^2)';    % Condizione iniziale
tspan = linspace(0,T,Nt);       % Tempi di output
% Chiamata al solver ODE
[t,Phi] = ode45(@ConvFun,tspan,Phi0);
% Grafica
for it = 1:Nt
    plot(x,Phi0,'.-r',x,Phi(it,:),'.-k'); axis([0 L -0.5 1.5]); drawnow;
end




