clc; close all; clear;

%%

z=10000;
epsilonD=0.98;
etaC=0.88;
etaT=0.94;
Qf=43000000;
etaB=0.96;
Au=0.126;
etaADN=0.96;
cp=1005;
gamma=1.4;
R=287;

nbeta=100;
beta=linspace(1,50,nbeta);
betaT=beta.^((gamma-1)/gamma);

T4=[1200 1400 1600];
M=[0.6 0.9 1.2 1.5];


