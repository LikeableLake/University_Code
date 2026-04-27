clc;
clear;
close all;

%% dati

%forze nodali
f2=1;
f3=1;
f4=1;
f5=1;

%f5 e f6 sono incognite dei vincoli
u1=0;
u6=0;

%nei vincoli 1 e 6 ho incastri: gli spostamenti sono nulli
k1= 10;
k2=20;
k3=100;
k4=20;
k5=10;

%% ricavo le matrici locali e globale

K1=zeros(6);
K2=zeros(6);
K3=zeros(6);
K4=zeros(6);
K5=zeros(6);

K1(1,:)=[k1 -k1 0 0 0 0];
K1(2,:)=[-k1 k1 0 0 0 0];

K2(2, :)=[0 k2 -k2 0 0 0];
K2(3, :)=[0 -k2 k2 0 0 0];

K3(3, :)=[0 0 k3 -k3 0 0];
K3(4, :)=[0 0 -k3 k3 0 0];

K4(4, :)=[0 0 0 k4 -k4 0];
K4(5, :)=[0 0 0 -k4 k4 0];

K5(5, :)=[0 0 0 0 k5 -k5];
K5(6, :)=[0 0 0 0 -k5 k5];

Kglob=K1+K2+K3+K4+K5;

disp('matrice delle rigidezze globale=');
disp(Kglob);

%% equazioni ricerca degli spostamenti

%visto che gli spostamenti u1 e u6 sono nulli posso ricavare una matrice
%ausiliaria che risolva il sistema solo per u2, u3, u4, u5

% Kaus=zeros(4);
% Kaus(1, :)=Kglob(2, 2:5);
% Kaus(2, :)=Kglob(3, 2:5);
% Kaus(3, :)=Kglob(4, 2:5);
% Kaus(4, :)=Kglob(5, 2:5);

Kaus=Kglob(2:5 , 2:5);

Faus=[f2;
    f3;
    f4;
    f5];

Uaus=Kaus\Faus;

%% risultati finali

U=[u1 ; Uaus; u6];

disp('vettore degli spostamenti nodali(in metri)=');
disp(U);

%% in aggiunta trovo le rezioni in 1 e 6

Fris=Kglob*U;

disp('vettore delle forze nodali(in Newton)=');
disp(Fris);

