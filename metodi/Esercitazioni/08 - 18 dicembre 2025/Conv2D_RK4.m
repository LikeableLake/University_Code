clc; close all; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che simula la equazione di convezione lineare 2D con un
% campo di velocità assegnato.
% La equazione è:
% 
%                 phi_t + u*phi_x + v*phi_y = 0
% 
% Dove u e v sono le componenti del campo di velocità, assegnato come un 
% moto di rotazione rigido. 
% Il codice adotta un approccio semidiscretizzato, nel quale la equazione 
% viene prima discretizzata nello spazio, e il risultante sistema di ODE 
% viene poi integrato nel tempo con un solver dedicato. In questa versione 
% viene implementato una classica procedura RK4.

global Nx Ny U V h         % Dichiarazione globale delle variabili

Lx = 1;  Ly = 1;           % Estensione del dominio
Nx = 80; Ny = 80;          % Numero di punti di discretizzazione
x  = linspace(-Lx,Lx,Nx);  % Mesh nella direzione x
y  = linspace(-Ly,Ly,Ny);  % Mesh nella direzione y
hx = x(2)-x(1);            % Passo spaziale lungo x
hy = y(2)-y(1);            % Passo spaziale lungo x
h  = hx;                   % Assumiamo passo spaziale uniforme e costante
% Campo di velocità ( rotazione rigida)
w  = 1;                    % Velocità angolare
u  = w*y;                  % Componente u
v  = w*(-x)';              % Componente v
U  = repmat(u,Nx,1);       % Matrice delle componenti u
V  = repmat(v,1,Ny);       % Matrice delle componenti u
% 
% Si osservi che gli array 2D contenenti le componenti di velocità (e più
% in generale tutti gli array 2D contenenti una funzione scalare) sono
% definiti secondo la convenzione per la quale l'indice i varia lungo la
% coordinata x e l'indice j lungo la coordinata y. 
% Questa convenzione è illustrata dallo schema:
%
%            ^
%            |
%         y  |                 f(i,j+1)
%            |                    o
%            |                    |
%            |                    |
%            |    f(i-1,j) o------o------o f(i+1,j)
%            |                  f(i,j)
%            |                    |
%            |                    o
%            |                 f(i,j-1)
%            |
%            |--------------------------------------------> 
%                                                       x
% 
% ed è diversa da quella adottata nelle routines grafiche di MATLAB, nelle
% quali il primo indice (i) varia lungo la direzione verticale (dal basso
% all'alto), mentre il secondo indice (j) varia lungo la direzione
% orizzontale (da sinistra a destra). 
% Per disegnare correttamente gli array 2D usando la grafica di MATLAB, è
% necessario trasporre gli array.
% 
% Grafica del campo di velocità
figure; quiver(x,y,U',V'); axis image; title('Campo di velocita'''); 
% Calcoliamo il Dt sulle max componenti di velocità
umax = max(max(abs(U))); vmax = max(max(abs(V)));
Cx   = 1.5;              Cy   = 1.5;
Dt   =  max([hx*Cx/umax,hy*Cy/vmax]);
T    =  2*pi/abs(w);           % Tempo finale di un 'giro'
% Modifichiamo il Dt e il Courant locale in modo che il tempo finale è 
% raggiunto dopo un numero intero di passi temporali 
Nt   =  round(T/Dt);     Dt = T/Nt;
% Condizioni iniziali
x0   = -Lx/2;   y0 = 0;        % Coordinate dell'origine della Gaussiana.
sig  = Lx/5;    sigq = sig^2;  % Varianza della Gaussiana.
PHI  = zeros(Ny,Nx);
for i = 1:Nx
    for j = 1:Ny
        PHI(i,j) = exp(-((x(i)-x0)^2 + (y(j)-y0)^2)/sigq);
    end
end
% Plot della condizione iniziale
figure(2); pcolor(PHI'); axis equal; shading interp; hold on; 
contour(PHI','k');  title('Condizione iniziale'); drawnow;  hold off
figure(3); surfl(PHI'); shading interp; hold on; 
contour3(PHI','k'); title('Condizione iniziale'); drawnow;  hold off
ax = axis;  pause;
% Ciclo instazionario
PHI_new = zeros(Ny,Nx);      % Allocazione
for it = 1:Nt
    t = it*Dt;
% Valori intermedi del RK
    PHI1 = PHI;              F1 = Conv2DFun(t,PHI1);
    PHI2 = PHI + 0.5*Dt*F1;  F2 = Conv2DFun(t,PHI2);
    PHI3 = PHI + 0.5*Dt*F2;  F3 = Conv2DFun(t,PHI3);
    PHI4 = PHI + Dt*F3;      F4 = Conv2DFun(t,PHI4);
% Step finale del RK
    PHI  = PHI + Dt*((1/6)*(F1+F4)+(1/3)*(F2+F3));
% Grafica
    figure(2); pcolor(x,y,PHI'); axis equal; shading interp
    hold on; contour(x,y,PHI','k'); hold off
    title(['t = ',num2str(t)]);
    drawnow;  
    figure(3); surfl(x,y,PHI'); 
    shading interp; hold on; contour3(x,y,PHI','k'); hold off
    axis([-Lx Lx -Ly Ly -0.2 1]); title(['t = ',num2str(t)]);
    drawnow
end











