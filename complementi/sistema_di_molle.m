clc; clear; close all;

%% ricerca dati
disp('calcolatore di spostamenti di sistemi di molle in serie')
N=input('inserisci il numero di molle \n');

while (isempty(N)||~isscalar(N)||floor(N)~=N||N<0)
    N=input('sbagliato, inserisci di nuovo \n');
end

%%
CostantiMolle=zeros(N, 1);
%%

for i=1:N
    CostantiMolle(i , 1)=input(sprintf('inserisci la costante elastica della molla numero %d (in N/m) \n', i));
    while (isempty(CostantiMolle(i, 1))||~isscalar(CostantiMolle(i, 1))||CostantiMolle(i,1)<0)
        CostantiMolle(i, 1)=input('sbagliato, inserisci di nuovo \n');
    end
end
%%
Fnodi=zeros(N-1, 1);
for i=1:N-1
    Fnodi(i)=input(sprintf('inserire la forza applicata al nodo %d in N (SOLO NODI INTERNI) \n', i));
    while (isempty(Fnodi(i))||~isscalar(Fnodi(i)))
        Fnodi(i)=input('sbagliato, inserisci di nuovo \n');
    end
end

%%  riepilogo dati e costruzione matrice di rigidezza

disp('-----------riepilogo forze nodali-----------')
disp('(le forze negative sono rivolte verso sinistra)')
disp(' ')

for i=1:N-1
    fprintf('f%d=%0.3f N \n', i, Fnodi(i))
end

disp(' ')

Kglob=zeros(N+1);

for i=1:N
    K=zeros(N+1);
    K(i:i+1, i:i+1)=[CostantiMolle(i) -CostantiMolle(i);
        -CostantiMolle(i) CostantiMolle(i)];
    Kglob=Kglob+K;
end

disp('la matrice di rigidezza è:')
disp(Kglob)

%% risoluzione problema

Kaus=Kglob(2:N, 2:N);
Uint=Kaus\Fnodi;

Utot=[0; Uint; 0];

disp(' ')
disp('ecco gli spostamenti di tutti i nodi(in m)')
disp(Utot)


%% reazioni vincolarri

Rvinc(1)=Kglob(1, :)*Utot;
Rvinc(2)=Kglob(N+1, :)*Utot;


disp('ecco le reazioni vincolari (in N)')
disp(Rvinc')
