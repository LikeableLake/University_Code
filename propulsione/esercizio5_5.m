clc; close all; clear;

%%

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
V=M.*a;

%% inizializzazione matrici delle soluzioni
S=zeros(nT, nbeta,nM);
TSFC=zeros(nT, nbeta,nM);
I=zeros(nT, nbeta, nM);
etaTH=zeros(nT, nbeta, nM);
etaP=zeros(nT,nbeta, nM);
eta=zeros(nT,nbeta, nM);

%% ciclo
for k=1:nM
    for j=1:nT
        for i=1:nbeta
            %presa d'aria
            if M(k)<=1
                eD=eDsubsonico;
            else
                eD=1-0.075*(M(k)-1)^1.35;
            end
            T2=Ta*(1+((gamma-1)/2)*M(k)^2);
            P2=eD*Pa;
            %compressore
            P3=beta(i)*P2;
            T3=T2*(1+(betaT(i)-1)/etaC);
            %combustione
            P4=P3;
            f=(cp*(T4(j)-T3))/(etaB*Qf);
            %turbina
            T5=T4(j)-(T3-T2)/(1+f);
            T5ideale=T4(j)-(T4(j)-T5)/etaT;
            P5=P4*((T5ideale/T4(j))^((gamma)/gamma-1));
            %ugello
            pressratio=Pa/P5;
            if (pressratio<=0.5283)
                %ugello strozzato
                Pu=P5*0.5283;
                Tu=T5*0.833;
                RHOu=Pu/(R*Tu);
                Vu=sqrt(gamma*R*Tu);
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            else 
                %l'ugello non è strozzato
                Pu=P5*pressratio;
                Tu=T5*((pressratio)^((gamma-1)/(gamma)));
                RHOu=Pu/(R*Tu);
                Vu=sqrt(2*cp*T5*etaADN*(1-(pressratio^((gamma-1)/gamma))));
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            end
            
            %risultati finali
            
            S(j,i,k)=Ma*((1+f)*Vu-V(k))+(Pu-Pa)*Au;
            I(j,i,k)=S(j,i,k)/Ma;
            TSFC(j,i,k)=f*Ma/S(j,i,k);
            etaTH(j,i,k)=(Vu^2-V(k)^2)/(2*f*Qf);
            etaP(j,i,k)=2*((V(k)/Vu)/(1+(V(k)/Vu)));
            eta(j,i,k)=etaTH(j,i,k)*etaP(j,i,k);

            %controllo se ci sono risultati immaginari, in questo caso il motore
            %non produce spinta infatti il rapporto di pressioni Pa/P05 è maggiore
            %di 1
            if(imag(Vu)~=0 || (real(S(j,i,k))<=0))
                S(j,i,k)=0;
                I(j,i,k)=0;
                TSFC(j,i,k)=inf;
                etaTH(j,i,k)=NaN;
                etaP(j,i,k)=NaN;
                eta(j,i,k)=NaN;
            end
        end
    end
end

    o=1;

for p=1:nM
    figure(o)
    plot(beta, I(1,:,p)); hold on;
    plot(beta, I(2,:,p)); hold on;
    plot(beta, I(3,:,p))
    title('Spinta specifica')
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    figure(o+1)
    plot(beta, TSFC(1,:,p)); hold on;
    plot(beta, TSFC(2,:,p)); hold on;
    plot(beta, TSFC(3,:,p))
    title('TSFC')
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    figure(o+2)
    plot(beta, etaTH(1,:,p)); hold on;
    plot(beta, etaTH(2,:,p)); hold on;
    plot(beta, etaTH(3,:,p))
    title('\eta_T_H')
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    figure(o+3)
    plot(beta, etaP(1,:,p)); hold on;
    plot(beta, etaP(2,:,p)); hold on;
    plot(beta, etaP(3,:,p))
    title('\eta_P')
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    figure(o+4)
    plot(beta, eta(1,:,p)); hold on;
    plot(beta, eta(2,:,p)); hold on;
    plot(beta, eta(3,:,p))
    title('\eta')
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    o=o+5;
end
