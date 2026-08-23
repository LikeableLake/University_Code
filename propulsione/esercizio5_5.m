clc; close all; clear;

%% svolgimento dell'esercizio 5.5
% l'esercizio è stato svolto utilizzando i valori dell'atmosfera isa per il
% calcolo di temperatura e pressione, la formula empirica per l'efficienza
% della presa d'aria supersonica, e un triplo ciclo for innestato in modo
% da ricavare le prestazioni in funzione delle tre variabili beta, T4, M0.

% L'utilizzo di matrici tridimensionali è stato utile al fine di conservare
% tutti i valori e di poter inserire in dei grafici le curve delle
% prestazioni in funzione di beta e T4 per poi analizzare il funzionamento
% in 4 mach diversi, uno per ogni figura che restituirà il codice.

z=10000;
eDsubsonico=0.98;
etaC=0.88;
etaT=0.94;
Qf=43000000;
etaB=0.96;
Au=0.126;
etaADN=0.96;
cp=1005;
gamma=1.4;
R=287;

nbeta=50;
beta=linspace(1,50,nbeta);
betaT=beta.^((gamma-1)/gamma);

T4=[1200 1400 1600];
nT=length(T4);

M=[0.6 0.9 1.2 1.5];
nM=length(M);


%% calcolo dei valori atmosferici
Ta=288-0.0065*z;
Pa=101330*((Ta/288)^(9.81/(R*0.0065)));
a=sqrt(gamma*R*Ta);
V0=M.*a;

%% inizializzazione matrici delle soluzioni
S=zeros(nT, nbeta,nM);
TSFC=zeros(nT, nbeta,nM);
I=zeros(nT, nbeta, nM);
etaTH=zeros(nT, nbeta, nM);
etaP=zeros(nT,nbeta, nM);
eta=zeros(nT,nbeta, nM);
adattato=ones(nT, nbeta, nM);

%% ciclo
for k=1:nM
    %presa d'aria
    if M(k)<=1
        eD=eDsubsonico;
    else
        eD=1-0.075*(M(k)-1)^1.35;
    end
    T2=Ta*(1+((gamma-1)/2)*M(k)^2);
    P0a = Pa*(1 + ((gamma-1)/2)*M(k)^2)^(gamma/(gamma-1));
    P2 = eD*P0a;

    for j=1:nT
        for i=1:nbeta
            %compressore
            P3=beta(i)*P2;
            T3=T2*(1+(betaT(i)-1)/etaC);
            %combustione
            P4=P3;
            f=(cp*(T4(j)-T3))/(etaB*Qf);
            %turbina
            T5=T4(j)-(T3-T2)/(1+f);
            T5ideale=T4(j)-(T4(j)-T5)/etaT;
            P5=P4*((T5ideale/T4(j))^((gamma)/(gamma-1)));
            %ugello
            pressratio=Pa/P5;
            if (pressratio<=0.5283)
                %ugello strozzato non adattato
                adattato(j,i,k)=0;
                Pu=P5*0.5283;
                Tu=T5*0.833;
                RHOu=Pu/(R*Tu);
                Vu=sqrt(gamma*R*Tu);
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            else                         
                %l'ugello non è strozzato quindi è adattato
                adattato(j,i,k)=1;
                Pu=P5*pressratio;
                Tu=T5*((pressratio)^((gamma-1)/(gamma)));
                RHOu=Pu/(R*Tu);
                Vu=sqrt(2*cp*T5*etaADN*(1-(pressratio^((gamma-1)/gamma))));
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            end
                   
            %risultati finali
            %utilizzo le formule generali in modo che anche con ugelli non
            %adattati i risultati siano coerenti            
                
            %calcolo la potenza del getto con la formula generale, con
            %effetto della propulsione ed effetto della velocità

            PJ = 0.5*Ma*((1+f)*Vu^2 - V0(k)^2)+ (Pu-Pa)*Au*Vu;

            %faccio la stessa cosa per la spinta e poi trovo la spinta
            %specifica e il TSFC

            S(j,i,k)=Ma*((1+f)*Vu-V0(k))+(Pu-Pa)*Au;
            I(j,i,k)=S(j,i,k)/Ma;
            TSFC(j,i,k)=f*Ma/S(j,i,k)*3600;
            
            %calcolo dei rendimenti

            etaP(j,i,k)=(S(j,i,k)*V0(k))/PJ;
            etaTH(j,i,k)=PJ/(Ma*f*Qf);
            eta(j,i,k)=etaTH(j,i,k)*etaP(j,i,k);

            %correzioni per situazioni fisicamente impossibili

            if Vu <= V0(k) || f<=0 || imag(Vu)~=0 || pressratio >= 1
                S(j,i,k)=NaN;
                I(j,i,k)=NaN;
                TSFC(j,i,k)=NaN;
                etaTH(j,i,k)=NaN;
                etaP(j,i,k)=NaN;
                eta(j,i,k)=NaN;
            end
        end
    end
end

%% stampa a video dei risultati

for p=1:nM

    figure(p)

    subplot(3, 2, 1)
    plot(beta, S(:,:,p)); hold on;
    title(sprintf('Spinta (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, 2)
    plot(beta, I(:,:,p)); hold on;
    title(sprintf('Spinta specifica (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, 3)
    plot(beta, TSFC(:,:,p)); hold on;
    title(sprintf('TSFC (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, 4)
    plot(beta, etaTH(:,:,p)); hold on;
    title(sprintf('\x03B7_T_H (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, 5)
    plot(beta, etaP(:,:,p)); hold on;
    title(sprintf('\x03B7_P (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, 6)
    plot(beta, eta(:,:,p)); hold on;
    title(sprintf('\x03B7 (M=%.1f)', M(p)))
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');
end
