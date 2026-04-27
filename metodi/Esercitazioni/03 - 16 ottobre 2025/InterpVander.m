clc; clear all; close all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Script di esercitazione sulla costruzione e sull'utilizzo della matrice
% di Vandermonde per la interpolazione Lagrangiana.

% Il presente codice effettua la interpolazione Lagrangiana 1D di dati 
% assegnati su mesh arbitrario costruendo la matrice di Vandermonde e
% risolvendo il sistema lineare.

% Assegnazioni preliminari.
g  = 5;                 % Grado del polinomio interpolante.
Np = g+1;               % Numero di punti del mesh.
L  = 2*pi;              % Definizione della lunghezza del dominio.
x  = linspace(0,L,Np);  % Mesh equispaziato sul dominio [0.L].
x  = x(:);              % Riduzione 'universale' a vettore colonna.
V  = nan(Np);           % Allocazione della matrice di Vandermonde.
f  = x.*sin(x);         % Valori della funzione di esempio da interpolare.

% Costruzione della matrice di Vandermonde.
for i = 1:Np
   V(:,i)= x.^(i-1);    % Assegnazione alla colonna i-ma di V.
end

% Calcolo dei coefficienti della interpolante
a  = V\f;               % Risoluzione del sistema lineare.
xf = linspace(0,L,100); % Mesh 'fitto' per la grafica.
ff = xf.*sin(xf);       % Funzione da interpolare sul mesh 'fitto'.
% Calcolo del polinomio interpolante sul mesh 'fitto'.
P = 0;
for j = 1:Np
    P = P + a(j)*xf.^(j-1);
end
err = norm(P-ff,inf);
% Grafica
plot(x,f,'k.','markersize',16); hold on;
plot(xf,P,'r-',xf,ff,'b--')
title (['Interpolazione di Lagrange di grado ',num2str(g)]);
text (0.5,-4,['Errore massimo assoluto: ',num2str(err)]);

% title(['Errore massimo assoluto: ',num2str(err)])
legend('Dati da interpolare',['Interpolante di grado ',num2str(g)],...
       'Funzione analitica')
xlabel('x'); ylabel('f(x)'); axis([0 L -6 3])
drawnow

