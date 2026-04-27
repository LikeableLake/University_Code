clc; clear all; close all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Script di esercitazione per la visualizzazione dei polinomi di Lagrange
% su mesh assegnato.
% Il codice utilizza le functions 'polyval.m' e 'polyfit.m' predefinite in
% Matlab per il calcolo delle interpolanti.

% Assegnazioni preliminari.
g  = 5;                 % Grado dei polinomi di Lagrange.
Np = g+1;               % Numero di punti del mesh.
L  = 1;                 % Definizione della lunghezza del dominio.
x  = linspace(0,L,Np);  % Mesh equispaziato sul dominio [0.L].
xf = linspace(0,L,100); % Mesh 'fitto' per la grafica.
I = eye(Np);            % Costruzione della matrice identica.

% Costruzione e grafica dei polinomi di Lagrange.
for i = 1:Np
% Il polinomio i-mo di Lagrange viene calcolato mediante il comando  
% polyval(polyfit(...),...) in corrispondenza dei dati assegnati mediante 
% la i-ma colonna della matrice identica, corrispondente a un vettore di
% tutti zero tranne la i-ma componente che e' pari a 1.
   P = polyval(polyfit(x',I(:,i),g),xf);
% Grafica
   plot(xf,P,'-',x,I(:,i)','.r','markersize',16); hold on;
   plot(x,0*x,'ko-');
   xlabel('x'); ylabel('L(x)'); axis([0,1,-2,2])
   title(['Polinomio di Lagrange n. ',num2str(i)])
   hold off;
   pause;
end