clc; close all; clear;

%% esercizio sui motori turbogetto

esercizio=1;
switch esercizio
    case 1
        z=10000;
        M=0.8;
        T4=1500;
        Qf=43000000;
        cp=1.005;
        gamma=1.4;
        R=cp*1000*(gamma-1)/gamma;
        eD=0.95;
        etaC=0.88;
        etaPB=1;
        etaB=0.95;
        etaT=0.94;

        %fan
        etaF=etaC;
        betaF=2;
        BPR=5;
        betaCF=9;

        %jet
        betaCJ=18;
        betaTJ=betaCJ^((gamma-1)/gamma);

        %dati atmosfera isa
        Ta=288-0.0065*z;
        Pa=101330*((Ta/288)^(9.81/(R*0.0065)));
        a=sqrt(gamma*R*Ta);
        V=M*a;

end

%soluzione

T2=Ta*(1+((gamma-1)/2)*M^2);
P02 = Pa*(1 + ((gamma-1)/2)*M^2)^(gamma/(gamma-1));
P2 = eD*P02;

% caso turbojet

%compressore
P3=betaCJ*P2;
T3=T2*(1+(betaTJ-1)/etaC);

%combustione
P4=P3*etaPB;
fJ=(cp*(T4-T3))/(etaB*Qf);

%turbina
T5=T4-(T3-T2)/(1+fJ);
T5ideale=T4-(T4-T5)/etaT;
P5=P4*((T5ideale/T4)^((gamma)/(gamma-1)));

%ugello
PuJ=Pa;

RHOuJ=PuJ/(R*TuJ);
VuJ=sqrt(2*cp*T5*(1-((PuJ/P5)^((gamma-1)/gamma))));
MuJ=RHOuJ*Au*VuJ;
MaJ=MuJ/(1+fJ);

%soluzioni finali

PjJ = 0.5*MaJ*((1+fJ)*VuJ^2 - V^2);

SJ=MaJ*((1+fJ)*VuJ-V)+(PuJ-Pa)*Au;
IJ=SJ/MaJ;
TSFCJ=fJ*MaJ/SJ*3600;

etaPJ=(SJ*V)/PjJ;
etaTHJ=PjJ/(MaJ*fJ*Qf);
etaJ=etaTHJ*etaPJ;