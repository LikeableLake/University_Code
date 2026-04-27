clc
clear
close all

%% primo quesito

disp('----------Parte_1----------');

%dati

u1=0;
u4=0;

f2=-10;
f3=10;

k1=30;
k2=30;
k3=15;
k4=30;

% svolgimento

k23=k2+k3;      %trovo la rigidezza della molla equivalente alle due in parallelo

K1=zeros(4);
K23=zeros(4);
K4=zeros(4);    %inizializzo le matrici di rigidezza locali già in forma 4*4

K1(1:2, 1:2)=[k1 -k1; -k1 k1];
K23(2:3, 2:3)=[k23 -k23; -k23 k23];
K4(3:4, 3:4)=[k4 -k4; -k4 k4];          %do i valori alle rigidezze locali tenendo conto della loro posizione nel sistema

Kglob=K1+K23+K4;        %Kglob è la somma delle matrici di rigidezza locale, dato che ho già tenuto conto della loro posizione in righe e colonne, nel passaggio precedente


%risolvo il problema trovando la matrice ridotta e calcolando gli
%spostamenti

Krid=Kglob(2:3, 2:3);
Frid=[f2 ; f3];

Urid=Krid\Frid;         %trovo u2 e u3 dal sistema ridotto

U=[u1; Urid; u4];         %ricavo il vettore U completo

disp('il vettore degli spostamenti è:')
fprintf('u%d=%0.3f m \n', 1, U(1));
fprintf('u%d=%0.3f m \n', 2, U(2));
fprintf('u%d=%0.3f m \n', 3, U(3));
fprintf('u%d=%0.3f m \n\n', 4, U(4));

%% secondo quesito
disp('----------Parte_2----------');

%i valori f2, k1, k4 sono gli stessi del passaggio precedente
%chiamerò u3_2 il valore di u3 in questa parte e così anche u2_2 e k23_2

u3_2=0.12;
u2_2=0.06;

Urid_2=[u2_2; u3_2];
% 
% %ricavo k23_2 dalla prima equazione
% k23_2=(f2-30*Urid_2(1))/(Urid_2(1)-Urid_2(2));
% %ricavo f3_2 dalla seconda
% f3_2=-k23_2*Urid_2(1)+(30+k23_2)*Urid_2(2);
% 

Krid_2=(f2*u2_2+f2*u3_2)/(Urid_2'*Urid_2)

disp('la rigidezza del parallelo è:')
fprintf('k23=%0.6f N/m \n\n', k23_2)

disp('la forza incognita è:')
fprintf('f3=%0.3f N/m \n\n', f3_2)

%% terzo quesito

disp('----------Parte_3----------');

%i dati sono gli stessi della parte 1 tranne f3


k1_3=k1;
k4_3=k4;
k23_3=k23;

Kglob_3=Kglob;
Krid_3=Krid;

f3_3=linspace(10, 30, 5);
f2_3=zeros(1, 5);
f2_3(:,:)=f2;
%ragiono in maniera matriciale, le f sono le forze applicate, non sono più
%scalari ma vettori riga, perchè a ogni colonna corrisponde un valore di f3_3
%diverso

Frid_3=[f2_3; f3_3];
U=zeros(2, 5);

for i=1:5
    U(:, i)=Krid_3\Frid_3(:, i);
    fprintf('per f3=%0.3f N ho: \n u2=%0.6f m \n u3=%0.6f m\n\n', f3_3(i), U(1,i), U(2,i))
end
%il ciclo for itera il calcolo degli spostamenti e li stampa nella command window, per ogni f3 diversa

%creo due grafici, uno per f3-u2 e uno per f3-u3 su un'unica figura
figure(1);
plot( U(1,:),f3_3, '-ok'); hold on
plot( U(2,:),f3_3, '-or'); hold on
grid on;
legend('u2', 'u3', 'Location','southeast');
ylabel('f3', 'Rotation',0, 'FontSize',14);
xlabel('u', 'Rotation',0, 'FontSize',14);

fprintf(['i grafici sono lineari, infatti già dalla formula iniziale possiamo\n' ...
    'capire che, essendo un sistema lineare, aggiungere una combinazione lineare\n' ...
    'del vettore F o della matrice K, semplicemente dà come risultato una\n' ...
    'combinazione lineare del vettore U\n']);