clc; clear; close all

%% inserimento dei dati

%il sistema è composto da 4 molle, 2 delle quali in parallelo, le altre due
%in serie con il sistma in parallelo
%ho quindi 4 nodi totali, 2 interni a cui sono applicate le forze:
Nnod=4;

f2=1;
f3=1;


%i nodi esterni sono vincolati quindi:

u1=0;
u4=0;

%costanti delle molle

k1=10;
k2=20;
k3=10;
k4=10;

%% svolgimento

%le due molle in parallelo possono essere sostituite da una molla
%equivalente di costante elastica:
k23=k2+k3;

%definisco le matrici di rigidezza

K1=zeros(Nnod);
K23=zeros(Nnod);
K4=zeros(Nnod);

K1(1:2,1:2)=[k1 -k1; -k1 k1];
K23(2:3,2:3)=[k23 -k23; -k23 k23];
K4(3:4,3:4)=[k4 -k4; -k4 k4];

Kglob=K1+K23+K4;

%dato che ho le condizioni al contorno elimino dalla Kglob le righe e le
%colonne corrispondenti agli spostamenti nulli

Kaus=Kglob(2:3,2:3);

Faus=[f2;f3];
Uaus=Kaus\Faus;
U=[0;Uaus;0];

disp("le forze applicate ai nodi interni sono:")
fprintf('f%d = %0.3f N \n', 2, f2)
fprintf('f%d = %0.3f N \n\n', 3, f3)

disp('la matrice di rigidezza globale del sistema in serie è:')
disp(Kglob)

disp('gli spostamenti dei nodi interni sono:')
fprintf('u%d = %0.3f m \n', 2, U(2))
fprintf('u%d = %0.3f m \n\n', 3, U(3))


%posso anche calcolare le forze incognite agenti sui vincoli:

f1=Kglob(1, :)*U;
f4=Kglob(4, :)*U;
F=[f1;f2;f3;f4];

disp('le reazioni vincolari sono:')
fprintf('f%d = %0.3f N \n', 1, f1)
fprintf('f%d = %0.3f N \n', 4, f4)


