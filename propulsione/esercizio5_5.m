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
adattato=logical(false(nT, nbeta, nM));

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
            P02 = Pa*(1 + ((gamma-1)/2)*M(k)^2)^(gamma/(gamma-1));
            P2 = eD*P02;
            %compressore
            P3=beta(i)*P2;
            T3=T2*(1+(betaT(i)-1)/etaC);
            %combustione
            P4=P3;
            f=(cp*(T4(j)-T3))/(etaB*Qf);
            %turbina
            T5=T4(j)-(T3-T2)/(etaT*(1+f));
            T5ideale=T4(j)-(T4(j)-T5)/etaT;
            P5=P4*((T5ideale/T4(j))^((gamma)/(gamma-1)));
            %ugello
            pressratio=Pa/P5;
            if (pressratio<=(2/(gamma+1))^(gamma/(gamma-1)))
                funzionamentomotore=logical(true);
                %ugello strozzato
                Pu=P5*0.5283;
                Tu=T5*0.833;
                RHOu=Pu/(R*Tu);
                Vu=sqrt(gamma*R*Tu);
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            elseif (pressratio>=1)
                %il motore non funziona
                funzionamentomotore=logical(false);
                S(j,i,k)=nan;
                I(j,i,k)=nan;
                TSFC(j,i,k)=nan;
                etaTH(j,i,k)=NaN;
                etaP(j,i,k)=NaN;
                eta(j,i,k)=NaN;
            else
                funzionamentomotore=logical(true);               
                %l'ugello non è strozzato
                adattato(j,i,k)=logical(true);

                Pu=P5*pressratio;
                Tu=T5*((pressratio)^((gamma-1)/(gamma)));
                RHOu=Pu/(R*Tu);
                Vu=sqrt(2*cp*T5*etaADN*(1-(pressratio^((gamma-1)/gamma))));
                Mu=RHOu*Au*Vu;
                Ma=Mu/(1+f);
            end
            
            if funzionamentomotore
                if P5 <= Pa || Vu <= V(k) || f<=0
                    funzionamentomotore = false;

                    S(j,i,k)=NaN;
                    I(j,i,k)=NaN;
                    TSFC(j,i,k)=NaN;
                    etaTH(j,i,k)=NaN;
                    etaP(j,i,k)=NaN;
                    eta(j,i,k)=NaN;
                end
            end

            %risultati finali
            
            if(funzionamentomotore)
                S(j,i,k)=Ma*((1+f)*Vu-V(k))+(Pu-Pa)*Au;
                I(j,i,k)=S(j,i,k)/Ma;
                TSFC(j,i,k)=f*Ma/S(j,i,k)*3600;
                Pj = 0.5*Ma*((1+f)*Vu^2 - V(k)^2)+ (Pu-Pa)*Au*Vu;
                etaP(j,i,k)=(S(j,i,k)*V(k))/Pj;
                etaTH(j,i,k)=Pj/(Ma*f*Qf);
                eta(j,i,k)=etaTH(j,i,k)*etaP(j,i,k);
            end
        end
    end
end

    o=1;

for p=1:nM

    figure(p)

    subplot(3, 2, o)
    plot(beta, I(1,:,p)); hold on;
    plot(beta, I(2,:,p)); hold on;
    plot(beta, I(3,:,p))
    title(['Spinta specifica (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, o+1)
    plot(beta, TSFC(1,:,p)); hold on;
    plot(beta, TSFC(2,:,p)); hold on;
    plot(beta, TSFC(3,:,p))
    title(['TSFC (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, o+2)
    plot(beta, etaTH(1,:,p)); hold on;
    plot(beta, etaTH(2,:,p)); hold on;
    plot(beta, etaTH(3,:,p))
    title(['\eta_T_H (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, o+3)
    plot(beta, etaP(1,:,p)); hold on;
    plot(beta, etaP(2,:,p)); hold on;
    plot(beta, etaP(3,:,p))
    title(['\eta_P (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, o+4)
    plot(beta, eta(1,:,p)); hold on;
    plot(beta, eta(2,:,p)); hold on;
    plot(beta, eta(3,:,p))
    title(['\eta (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');

    subplot(3, 2, o+5)
    plot(beta, adattato(1,:,p),'+'); hold on;
    plot(beta, adattato(2,:,p),'o'); hold on;
    plot(beta, adattato(3,:,p),'*');
    title(['l''ugello è adattato? (M=' num2str(M(p)) ')'])
    legend('T4=1200 K', 'T4=1400 K', 'T4=1600 K');
    ylim([-0.1 1.1]);
    text(40, 0.5, sprintf('0 = NO\n1 = SI'))

end
