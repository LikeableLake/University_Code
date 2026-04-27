clear all; clc; close all

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
%
% Script di test per la function PesiDer.
% Calcola la derivata di una funzione in un punto ed effettua uno studio di 
% accuratezza delle formule.

% Definiamo un vettore hv di passi spaziali per il calcolo dell'errore delle
% formule al variare di h.
H  = logspace(-5,1,10);
N  = length(H);                    % Numero di passi h da valutare.
v  = -1:1;                         % Indici dello stencil.
d  = 2;                            % Grado di derivazione richiesto
for i = 1:N
    h   = H(i);                    % Passo spaziale
    xs  = v*h;                     % Stencil
    xc  = xs(2);                   % Punto di collocazione
    w   = PesiDer(xs,xc,d);        % Calcolo i pesi.
    f   = exp(xs);                 % Vettore dei termini noti.
    dfn = w*f';                    % Derivata numerica.
    df  = exp(xc);                 % Derivata analitica.
    err(i) = abs(df-dfn);          % Calcolo dell'errore della formula.
end
% Grafica
loglog(H,err,'k.-','markersize',12);  grid on; hold on;
loglog(H,H.^2,'r-');
xlabel('\Delta x');   ylabel('Errore')

    


