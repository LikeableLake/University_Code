close all; clear all; clc;

%% SOLUZIONE DELLA TRACCIA 20/02/2013 - A

% 1) Risolvere l'equazione di laplace ∇²ζ = 0 su dominio quadrato con BC
% alla Dirichlet,omogenee su 3 lati e non omogenea sul lato restante

% 2) Risolvere l'equazione di Poisson ∇²ψ = ζ su dominio quadrato con
% BC alla Dirichlet omogenee sui 4 lati

% 3) Risolvere il punto 1) su dominio rettangolare, assegnando due
% BC non omogenee per ζ sui lati "corti" del rettangolo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4) Esercizio aggiuntivo: Risolvere il punto 3) Con BC alla Dirichlet
% omogenee sui lati "lunghi", una BC alla Dirichlet non omogenea su un lato
% corto e una Neumann omogenea sul lato restante. Usiamo ζ del punto 3)
% come termine sorgente

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTA: Nel codice seguente è spesso non necessario ridefinire alcune 
% grandezze fissate nei punti precedenti; ridefiniamo tutto per ogni punto
% della traccia per chiarezza


%% PUNTO 1

N    = 640;                    % Numero di punti lungo una direzione
L    = 1;                     % Lunghezza del lato
x    = linspace(0,L,N);       % Discretizzazione delle direzioni coordinate
y    = x;
y    = fliplr(y);             % Orientiamo y verso l'alto
h    = x(2) - x(1);           % Mesh uniforme
hq   = h*h;
G    = numgrid('S',N);        % Generiamo la griglia ordinata
L    = -delsq(G)/hq;          % Operatore Laplaciano su mesh uniforme
Ninc = max(max(G));           % Numero di incognite
q    = zeros(N);              % Inizializzazione dei termini noti
ZetaB= 1;                     % Zeta sul bordo per la BC non omogenea

for i = 1:N
    for j = 2:N               % Parte da 2 per calcolare j-1 anche con j==0
        if ( G(i,j)~=0 && G(i,j-1)==0 )   % Colonna Ovest
            q(i,j) = q(i,j) - ZetaB/hq;   % Modifica dei termini noti
        end
    end
end

q = q(2:N-1,2:N-1);                       % Risolviamo per i punti interni
q = q(:);                                 % Mega-vettore (N-2)x(N-2)
Zeta = L\q;                               % Risoluzione del sistema lineare
Zeta = reshape(Zeta,N-2,N-2);             % Riversiamo nella matrice 
Zeta = [ZetaB*ones(N-2,1), Zeta, zeros(N-2,1)];% Orliamo con le BC
Zeta = [zeros(1,N); Zeta; zeros(1,N)];


%  Grafica
figure(1)
surf(x,y,Zeta);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\zeta$','Fontsize',16,'Interpreter','latex');
title('Soluzione di $\nabla^2\zeta=0$ con BC alla Dirichlet sul quadrato',...
    'Interpreter','latex');
axis('square');

figure(2)
contour(x,y,Zeta);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\zeta$','Fontsize',16,'Interpreter','latex');
title('Isolinee di $\nabla^2\zeta=0$ con BC alla Dirichlet sul quadrato',...
    'Interpreter','latex');
axis('square');

%% PUNTO 2

% Il termine di produzione sarà la soluzione del punto precedente
p = Zeta(2:N-1,2:N-1); 
p = p(:);
Psi = L\p;
Psi = reshape(Psi,N-2,N-2);
Psi = [zeros(N-2,1), Psi, zeros(N-2,1)];
Psi = [zeros(1,N); Psi; zeros(1,N)];


figure(3)
surf(x,y,Psi);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\zeta$','Fontsize',16,'Interpreter','latex');
title('Soluzione di $\nabla^2\psi=\zeta$ con BC alla Dirichlet sul quadrato',...
    'Interpreter','latex');
axis('square');

figure(4)
contour(x,y,Psi);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\psi$','Fontsize',16,'Interpreter','latex');
title('Isolinee di $\nabla^2\psi=\zeta$ con BC alla Dirichlet sul quadrato',...
    'Interpreter','latex');
axis('square');


%% PUNTO 3
Lx     = 1;
Ly     = 2*Lx;                % Dominio rettangolare
Nx     = N;
Ny     = 1 + 2*(Nx - 1);      % Dalla condizione hx=hy 
x      = linspace(0,Lx,Nx); 
y      = linspace(0,Ly,Ny); 
y      = fliplr(y);
h      = x(2) - x(1);
hq     = h*h;
G      = 1:(Nx-2)*(Ny-2);      % Costruiamo la matrice G
G      = reshape(G,Ny-2,Nx-2); % riversiamo in una matrice e orliamo con 0
G      = [zeros(Ny-2,1),G,zeros(Ny-2,1)];
G      = [zeros(1,Nx); G ; zeros(1,Nx)];
L      = -delsq(G)/hq;

q      = zeros(Ny,Nx);
ZetaB  = 1;                    % BC alla Dirichlet sui lati "corti"             
ZetaB2 = 1;

for i = 1:Nx
    for j = 2:Ny-1             % Definito anche per j==0 e j==Ny
        if (G(j,i)~=0 && G(j+1,i)==0)
            q(j,i) = q(j,i) - ZetaB/hq;
        elseif (G(j,i)~=0 && G(j-1,i)==0)
            q(j,i) = q(j,i) - ZetaB2/hq;
        end
    end
end

q   = q(2:end-1,2:end-1);
q   = q(:);
Zeta = L\q;
Zeta = reshape(Zeta,Ny-2,Nx-2);
Zeta = [zeros(Ny-2,1), Zeta, zeros(Ny-2,1)];
Zeta = [ZetaB*ones(1,Nx); Zeta; ZetaB2*ones(1,Nx)];


figure(5)
surf(x,y,Zeta);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\zeta$','Fontsize',16,'Interpreter','latex');
title(['Soluzione di $\nabla^2\zeta=0$ con BC alla Dirichlet '...
    'sul rettangolo'],'Interpreter','latex');
axis('square');

figure(6)
contour(x,y,Zeta);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\zeta$','Fontsize',16,'Interpreter','latex');
title(['Isolinee di $\nabla^2\zeta=0$ con BC alla Dirichlet '...
    'sul rettangolo'],'Interpreter','latex');
axis('square');



%% PUNTO AGGIUNTIVO
Zeta_p = Zeta; % La soluzione del punto precedente diventa la sorgente

Lx     = 1;
Ly     = 2*Lx;
Nx     = N;
Ny     = 1 + 2*(Nx - 1);
x      = linspace(0,Lx,Nx); 
y      = linspace(0,Ly,Ny); 
y      = fliplr(y);
h      = x(2) - x(1);
hq     = h*h;
G      = 1:(Nx-2)*(Ny-2);
G      = reshape(G,Ny-2,Nx-2);
G      = [zeros(Ny-2,1),G,zeros(Ny-2,1)];
G      = [zeros(1,Nx); G ; zeros(1,Nx)];
L      = -delsq(G)/hq;

q      = Zeta_p;

PsiB  = 1;
g      = -1;        % Valore della derivata direzionale sul bordo (Neumann)

for i = 1:Nx
    for j = 2:Ny-1  % Definito anche per j==0 e j==Ny
        if (G(j,i) ~= 0 && G(j+1,i) == 0)
                q(j,i) = q(j,i) - PsiB/hq;
        end
        if (G(j,i) ~= 0 && G(j-1,i) == 0)
                node = G(j,i);
                L(node,node) = L(node,node) + 1/hq;
                q(j,i) = q(j,i) + g/h;    % Modifica dovuta alla BC Neumann
        end
    end
end


q    = q(2:end-1,2:end-1);
q    = q(:);
Psi = L\q;
Psi = reshape(Psi,Ny-2,Nx-2);
Psi = [zeros(Ny-2,1), Psi, zeros(Ny-2,1)];
Psi = [Psi(1,:); Psi; PsiB*ones(1,Nx)];



figure(7)
surf(x,y,Psi);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\Psi$','Fontsize',16,'Interpreter','latex');
title(['Soluzione di $\nabla^2\Psi=\zeta$ con BC alla Neumann/Dirichlet '...
    'sul rettangolo'],'Interpreter','latex');
axis('square');

figure(8)
contour(x,y,Psi);
colormap('hot'); colorbar
shading interp
xlabel('x','Fontsize',16,'Interpreter','latex');  
ylabel('y','Fontsize',16,'Interpreter','latex'); 
zlabel('$\Psi$','Fontsize',16,'Interpreter','latex');
title(['Isolinee di $\nabla^2\Psi=\zeta$ con BC alla Neumann/Dirichlet '...
    'sul rettangolo'],'Interpreter','latex');
axis('square');


