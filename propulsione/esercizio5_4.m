clc; close all; clear;

%%

Ta=288;
Pa=101330;

Ed=0.95;
etaC=0.88;
etaT=0.94;
Qf=43000000;
Au=0.126;
cp=1005;
gamma=1.4;
R=287;
T4=[1200 1400 1600];
nT=length(T4);

nbeta=50;

beta=linspace(1,50,nbeta);
betaT=beta.^((gamma-1)/gamma);

%% ciclo

S=zeros(nT, nbeta);
TSFC=zeros(nT, nbeta);
I=zeros(nT, nbeta);
etaTH=zeros(nT, nbeta);

%%
for j=1:nT
    for i=1:nbeta
    %% presa d'aria
    T02=Ta;
    P02=Ed*Pa;
    %% compressione
    P03=P02*beta(i);
    T3ideale=T02*betaT(i);
    T3reale=T02*(1+((betaT(i)-1)/etaC));
    %% combustione
    P04=P03;
    f=(cp*(T4(j)-T3reale))/Qf;
    %% turbina
    T5reale=T4(j)-(T3reale-T02)/(1+f);
    T5ideale=T4(j)-(T4(j)-T5reale)/etaT;
    P05=P04*((T5ideale/T4(j))^((gamma)/gamma-1));
    %% ugello
    pressratio=Pa/P05;
    if (pressratio<=0.5283)
        % l'ugello è strozzato
        Pu=P05*0.5283;
        Tu=T5reale*0.833;
        RHOu=Pu/(R*Tu);
        Vu=sqrt(gamma*R*Tu);
        Mu=RHOu*Au*Vu;
        Ma=Mu/(1+f);
    else 
        %l'ugello non è strozzato
        Pu=P05*pressratio;
        Tu=T5reale*((pressratio)^((gamma-1)/(gamma)));
        RHOu=Pu/(R*Tu);
        Vu=sqrt(2*cp*T5reale*(1-(pressratio^((gamma-1)/gamma))));
        Mu=RHOu*Au*Vu;
        Ma=Mu/(1+f);
    end
    %% risultati finali
    S(j,i)=Mu*Vu+(Pu-Pa)*Au;
    I(j,i)=S(j,i)/Ma;
    TSFC(j,i)=(Mu-Ma)/S(j,i);
    etaTH(j,i)=(Vu^2)/(2*f*Qf);

    %controllo se ci sono risultati immaginari, in questo caso il motore
    %non produce spinta infatti il rapporto di pressioni Pa/P05 è maggiore
    %di 1
    if(imag(Vu)~=0)
        S(j,i)=0;
        I(j,i)=0;
        TSFC(j,i)=inf;
        etaTH(j,i)=NaN;
    end
    end
end

figure (1);
plot(beta, I(1,:)); hold on;
plot(beta, I(2,:)); hold on;
plot(beta, I(3,:))
title('Spinta specifica')
legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

figure(2)
plot(beta, TSFC(1,:)); hold on;
plot(beta, TSFC(2,:)); hold on;
plot(beta, TSFC(3,:))
title('TSFC')
legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

figure(3)
plot(beta, etaTH(1,:)); hold on;
plot(beta, etaTH(2,:)); hold on;
plot(beta, etaTH(3,:))
title('\eta_T_H')
legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');