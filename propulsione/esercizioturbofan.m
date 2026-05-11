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

        %dati atmosfera isa
        Ta=288-0.0065*z;
        Pa=101330*((Ta/288)^(9.81/(R*0.0065)));
        a=sqrt(gamma*R*Ta);
        V=M.*a;

        %


end

