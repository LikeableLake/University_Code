clc;clear;close all;

% discretizzazione del dominio

Nx=50;
L=1;
x=linspace(0,L,Nx);
h=x(2)-x(1);
hq=h^2;

% definizione del tempo

beta=0.45;
k=1;
dt=hq*beta/k;
T=(L^2)*beta/k/4;
Nt=round(T/dt);
dt=T/Nt;

% condizione iniziale 

PHI0=1+x/L+sin(x*2*pi/L)+3*sin(x*4*pi/L);

figure(1);
plot(x, PHI0, 'r.-');
title('Condizone iniziale');
drawnow;
hold on;


% condizioni al contorno e termine produzione

P=zeros(Nx-2,1)*dt;

P(1)=P(1)+beta*1;
P(end)=P(end)+beta*2;

% P è costante nel tempo quindi non va discretizzata in n, n+1/2 e n+1

% creazione operatore derivata

D2=-gallery('tridiag', Nx-2);

% soluzione chiamo n+1/2 'ext'

PHI=PHI0(2:end-1)';

err=0.01;
DPHI=inf;
t=0;

%%

while DPHI/dt>err
    
    OLDPHI=PHI;

    PHIext=PHI+beta*D2*PHI/2+P/2;
    PHI=PHI+beta*D2*PHIext+P;
    
    DPHI=abs(max(-OLDPHI+PHI));

    plot(x, [1; PHI; 2], 'w:', x, PHI0, 'r.-');
    title(['soluzione all''equazione per t=' num2str(t)]);
    hold off;
    drawnow;
    t=t+dt;

end

%% seconda parte

I=eye(Nx-2);
Nb=20;
beta=linspace(0,1,Nb);
beta=beta(2:end);
error=zeros(Nb-1, 1);

for nbeta=1:Nb-1
    nbeta
    PHI=PHI0(2:end-1)';
    
    beta=0.45;
    k=1;
    dt=hq*beta/k;
    T=(L^2)*beta/k/4;
    Nt=round(T/dt);
    dt=T/Nt;

    err=0.01;
    DPHI=inf;
    t=0;

    while DPHI/dt>err
        OLDPHI=PHI;

        Trasf=I+beta(nbeta)*D2+((beta(nbeta))^2)*D2*D2/2;
        PHI=Trasf*PHI;

        DPHI=abs(max(-OLDPHI+PHI));
        t=t+dt;        

        error(nbeta)=norm(Trasf);
    end
end
%%
figure(2); plot(beta, error');
text(0.25, 1.2, 'l''instabilita'' parte da \beta = 0.5'); hold on;


%% terza parte
%riinizializzo i valori che sono cambiati nei punti precedenti

PHI=PHI0(2:end-1)';
err=0.01;
DPHI=inf;
t=0;
beta=0.45;
k=1;
dt=hq*beta/k;
T=(L^2)*beta/k/4;
Nt=round(T/dt);
dt=T/Nt;

Trasf=I+beta*D2+(beta^2)*D2*D2/2;

%soluzione analitica

PHIA=