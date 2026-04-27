%% 
clc; clear; close all;
%% disegno del problema
lenght=2;

nodes=[ 0       0;
        lenght/2     (((3)^(0.5))/2)*lenght; ...
        lenght       0; ...
        1.5*lenght   (((3)^(0.5))/2)*lenght];


figure(1); hold on; grid on;
xlabel('x'); ylabel('y');
axis('equal', 'padded');
plot(nodes(:,1), nodes(:,2), 'ow', 'MarkerFaceColor', 'w', 'MarkerSize', 10)
plot(nodes(:,1), nodes(:,2), '-g', 'LineWidth', 2)
plot([nodes(1,1);nodes(3,1)], [nodes(1,2),nodes(3,2)], '-g', 'LineWidth', 2)
plot([nodes(2,1);nodes(4,1)], [nodes(2,2),nodes(4,2)], '-g', 'LineWidth', 2)
plot([nodes(1,1);nodes(2,1)], [nodes(1,2),nodes(2,2)], 'squareb', 'MarkerFaceColor', 'b', 'MarkerSize', 20)
quiver(nodes(3,1), nodes(3,2), 0, lenght/2, 'r', 'LineWidth', 3, 'DisplayName', 'F');
quiver(nodes(4,1), nodes(4,2), lenght/2, 0, 'r', 'LineWidth', 3);


%% svolgimento esercizio
syms l E A real;
U=sym('u', [8, 1]);


%% asta 12

b12=sym (60/180)*pi;


c=cos(b12);
s=sin(b12);

K12loc=simplify([c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
       -c^2 -c*s c^2 c*s;
       -c*s -s^2 c*s s^2]);

F12loc=K12loc*U(1:4)*E*A/l;

%% asta 13

b13=sym (0/180)*pi;


c=cos(b13);
s=sin(b13);

K13loc=simplify([c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
       -c^2 -c*s c^2 c*s;
       -c*s -s^2 c*s s^2]);

F13loc=K13loc*U([1 2 5 6])*E*A/l;

%% asta 23

b23=sym (-60/180)*pi;


c=cos(b23);
s=sin(b23);

K23loc=simplify([c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
       -c^2 -c*s c^2 c*s;
       -c*s -s^2 c*s s^2]);

F23loc=K23loc*U(3:6)*E*A/l;

%% asta 24

b24=sym (0/180)*pi;


c=cos(b24);
s=sin(b24);

K24loc=simplify([c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
       -c^2 -c*s c^2 c*s;
       -c*s -s^2 c*s s^2]);

F24loc=K24loc*U([3 4 7 8])*E*A/l;



%% asta 34

b34=sym (60/180)*pi;


c=cos(b34);
s=sin(b34);

K34loc=simplify([c^2 c*s -c^2 -c*s;
        c*s s^2 -c*s -s^2;
       -c^2 -c*s c^2 c*s;
       -c*s -s^2 c*s s^2]);

F34loc=K34loc*U(5:8)*E*A/l;

%% sommo tutte le forze

F12=[F12loc; 0; 0; 0; 0];
F13=[F13loc(1:2); 0;0;F13loc(3:4); 0;0];
F23=[0;0;F23loc;0;0];
F24=[0;0; F24loc(1:2); 0;0; F24loc(3:4)];
F34=[0;0;0;0; F34loc];

Ftot=F12+F13+F23+F24+F34;

Ktot=simplify(jacobian(Ftot,U));

%% calcolo degli spostamenti dei nodi liberi

Fnote=sym([0; 1; 1; 0]);
Krid=Ktot(5:8,5:8);

Urid=simplify(inv(Krid)*Fnote);

U(5:8)=Urid;
U(1:4)=0;

%% calcolo delle reazioni vincolari

Fvinc=simplify(Ktot(1:4,5:8)*U(5:8));

%% esposizione finale

Ftot=[Fvinc; Fnote]
U