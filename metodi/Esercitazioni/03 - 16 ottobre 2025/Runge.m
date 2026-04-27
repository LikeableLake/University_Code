clc; clear all; close all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Script di esercitazione che implementa la interpolazione della funzione
% di Runge (f = 1/(1+x^2)) con i metodi:
% - Interpolazione Lagrangiana su mesh uniforme
% - Interpolazione SPLINE (cubica) su mesh uniforme.
% - Interpolazione Lagrangiana su mesh di Chebyshev.
% 

% Assegnazioni preliminari.
g  = 10;                  % Grado della interpolante Lagrangiana.
Np = g + 1;               % Numero di punti nel dominio.
xi = -4;     xf =  4;     % Estremi del dominio di interpolazione.
xg = linspace(xi,xf,200); % Mesh 'fitto' per la rappresentazione grafica.
fg = 1./(1+xg.^2);        % Funzione di Runge per la grafica.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolazione su mesh equispaziato (Lagrangiana e SPLINE).
x  = linspace(xi,xf,Np);  % Mesh di interpolazione equispaziato.
f  = 1./(1+x.^2);         % Dati della funz. di Runge su mesh equispaziato.
a  = polyfit(x,f,g);      % Coefficienti della interpolante su x.
fp = polyval(a,xg);       % Valori della interpolante su xg.
fs = spline(x,f,xg);      % Valori della SPLINE su xg.
dL = norm(fg-fp,inf);     % Errore massimo della interpolante Lagrangiana.
dS = norm(fg-fs,inf);     % Errore massimo della SPLINE.
% Grafica.
figure(1); subplot(3,1,1);
plot(xg,fg,'k-', xg,fp,'r-',x,f,'k.','markersize',16) 
title(['Interpolante di grado ',num2str(g),' su mesh equispaziato'])
text(-3.5,1.7,['Errore massimo: ',num2str(dL)])
xlabel('x'); ylabel('f(x)'); axis([xi,xf,-1,2]);
subplot(3,1,2);
plot(xg,fg,'k-', xg,fs,'r-',x,f,'k.','markersize',16) 
title(['Interpolante SPLINE cubica su ',num2str(g),' intervalli equispaziati'])
text(-3.5,1.7,['Errore massimo: ',num2str(dS)])
xlabel('x'); ylabel('f(x)'); axis([xi,xf,-1,2]);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interpolazione su mesh di Chebyshev
t   = linspace(0,pi,Np);  % Equispaziatura dell'angolo teta fra 0 e pi.
xch = xf*cos(t);          % Mesh di Chebyshev.
fch = 1./(1+xch.^2);      % Dati della funz. di Runge su mesh di Chebyshev.
a   = polyfit(xch,fch,g); % Coefficienti della interpolante su xch.
fp  = polyval(a,xg);      % Valori della interpolante su x.
d   = norm(fg-fp,inf);    % Errore massimo della interpolante.
% Grafica del confronto fra i mesh
figure(2) 
plot(x,0*x,'r.-','markersize',16); hold on;
plot(xch,0*xch+1,'k.-','markersize',16)
title(['Confronto fra mesh equispaziato e di Chebyshev su ',...
    num2str(Np),' punti']);  axis([-5 5 -1 2])
% Grafica
figure(1)
subplot(3,1,3);
plot(xg,fg,'k-', xg,fp,'r-',xch,fch,'k.','markersize',16)
title(['Interpolante di grado ',num2str(g),' su mesh di Chebyshev'])
text(-3.5,1.7,['Errore massimo: ',num2str(d)])
xlabel('x'); ylabel('f(x)'); axis([xi,xf,-1,2]);
