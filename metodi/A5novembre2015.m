clc; clear; close all;

% discretizzazione quadrato

N=40;
L=1;
x=linspace(0,L,N);

G=numgrid('S', N);
h=x(2)-x(1);
hq=h^2;

% definizione laplaciano

L=-delsq(G);

% definizione termine sorgente

P=ones(N-2);

% condizioni al contorno dirchlet = 0 ovunque quindi dovrei aggiungere 0 al
% contorno di P

%risoluzione

L=L/hq;

PHIPOISSON=L\P(:);

solp=zeros(N,N);
solp(2:N-1,2:N-1)=reshape(PHIPOISSON, N-2, N-2);

figure(1); subplot(1,2,1); contour3(x, x, solp, 'k'); hold on;
surf(x,x,solp, LineStyle= 'none', FaceColor='interp');

title('soluzione \phi all''equazione di poisson');
figure(1); subplot(1,2,2); contour(x,x,solp);  axis square;
title('grafico di \phi con le contour lines'); colorbar;

%%parte 2

%ridefinisco la condizione sorgente e il laplaciano

L=-delsq(G);
P=zeros(N-2);

% le condizioni al contorno sono tutte dirichlet nulle quindi non devo
% modificare ne' P ne' L

%definisco theta e beta

theta=0;
beta=1;

%definisco le matrici che moltiplicano phi(n+1) e phi(n)

I=eye(G(end-1,end-1));
B=I-beta*theta*L;
A=I+beta*(1-theta)*L;


%definisco il tempo

k=0.1;
T=0.25;
dt=beta*hq/k;
Nt=round(T/dt);
dt=T/Nt;

%condizione iniziale

PHIFOURIER=PHIPOISSON;

%% risolvo l'equazione


for i=1:Nt

    PHIFOURIER=B\(A*PHIFOURIER);

    solf=zeros(N,N);
    solf(2:N-1,2:N-1)=reshape(PHIFOURIER, N-2, N-2); 
    
    figure(2); surf(x,x,solf); shading interp; hold on; 
    contour3(x,x,solf,'k'); hold off; zlim([-0.1 0]);
    title('soluzione \phi all''equazione di fourier'); 
    text(0.5, 0.5, -0.08, ['tempo passato= ' num2str(i/Nt) ' s'] ); drawnow;

end


%% calcolo dell'errore



ntheta=10;
nbeta=10;
error=zeros(ntheta, nbeta);
theta=linspace(0,1,ntheta);
beta=linspace(0,1,nbeta);

for itheta=1:ntheta 
    itheta
    for ibeta=1:nbeta
       
        I=eye(G(end-1,end-1));
        B=I-beta(ibeta)*theta(itheta)*L;
        A=I+beta(ibeta)*(1-theta(itheta))*L;
       
        dt=beta(ibeta)*hq/k;
        Nt=round(T/dt);
        dt=T/Nt;

        Trasf=B\A;

        error(itheta,ibeta)= norm(Trasf);

    end
end


figure(3); subplot(1,2,1); contour3(beta, theta,  error, 'k'); hold on;
surf(theta,beta,error, LineStyle= 'none', FaceColor='interp');

title('norma della matrice di trasformazione, stabilità dello schema');
ylabel('\theta'); xlabel('\beta');

figure(3); subplot(1,2,2); contour(beta, theta, error);  axis square;
title('contour lines'); ylabel('\theta'); xlabel('\beta');
