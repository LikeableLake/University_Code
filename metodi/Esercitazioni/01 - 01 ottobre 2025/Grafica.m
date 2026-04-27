close all; clear all; clc;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
%
% Script di prova di comandi grafici in MATLAB


% Grafica instazionaria 1D
N  = 30;                      % Numero di punti
a  = 0;       b = 2*pi;       % Estremi dell'intervallo
x  = linspace(a,b,N);         % Mesh spaziale
T  = pi;                      % Tempo finale
Dt = 0.01;                    % Intervallo temporale
Nt = round(T/Dt);             % Numero di Dt
t  = 0;                       % Tempo iniziale
for i = 1:Nt
    y = sin(x-2*t); 
    t = t+Dt;
    figure(1); plot(x,y,'.-'); title('Onda che viaggia');
    xlabel('x'); ylabel('y');   drawnow;
end
