%%esercizio lez4 6/10/2025

clc;
clear;
close all;

%%

k1=100;
k2=100;

F=[0;10;0];

fprintf('f1=%.3f N\n', F(1));
fprintf('f2=%.3f N\n', F(2));
fprintf('f3=%.3f N\n\n', F(3));

K1=zeros(3);
K2=zeros(3);
U=zeros(3,1);

K1(1,:)=[k1 -k1 0];
K1(2,:)=[-k1 k1 0];

K2(2,:)=[0 k2 -k2];
K2(3,:)=[0 -k2 k2];

K=K1+K2;


disp ('matrice di rigidezza K=')
disp (K)

U(2)=F(2)/K(2,2);

fprintf('u1=%.3f m\n', U(1));
fprintf('u2=%.3f m\n', U(2));
fprintf('u3=%.3f m\n\n', U(3));


%% grafico f2 u2 con ciclo for



Ftot=zeros(1,10);
Utot=zeros(1,10);

Ftot(1)=F(2);
Utot(1)=U(2);

for i=1:9
    F(2)=F(2)+1;
    U(2)=F(2)/K(2,2);
    Ftot(i+1)=F(2);
    Utot(i+1)=U(2);
end

plot(Utot,Ftot);

%% grafico a mano e interpolando

X=[10; 11; 12; 13; 14; 15; 16; 17;18;19];
Xinterp=interp(X,10);

Y=Xinterp/K(2,2);

plot(Xinterp,Y)

%% grafico usando linspace


F2=linspace(0, 20, 21);
U2=F2/200;


U2interp=interp(U2, 2);  %nel dominio interpolo per trovare altri punti, l'interpolazione è lineare e usa i punti di riferimento F2 e U2
F2interp=interp(F2, 2);

figure(3);
plot(U2interp,F2interp, 'k.');
hold on;
plot(U2,F2, 'ro');
legend('interpolati', 'non interpolati');
grid on;
xlabel('spostamenti');
ylabel('forza');
title('molla');



