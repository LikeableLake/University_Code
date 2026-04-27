function TreMasse

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
%
% Applicazione degli integratori di ODE al problema dei tre corpi.
% Il codice scrive e risolve il sistema di equazioni differenziali
% ordinarie che governa il moto di tre masse ciascuna sotto l'influenza del
% campo gravitazionale delle altre due.
% 
clc; close all;

global m1 m2 m3

% Masse dei tre corpi.
m1 = 500;     m2 = 2;        m3 = 1;
% Posizioni iniziali dei tre corpi.
x1 = 0;        x2 = 2;        x3 = 0;        
y1 = 0;        y2 = 0;        y3 = 1;
% Velocita' iniziali dei tre corpi.
u1 = 0;        u2 = 0;        u3 = -20;        
v1 = 0;        v2 = 10;       v3 = 0;
% Assemblamento vettore delle condizioni iniziali.
y0 = [x1,y1, u1,v1, x2,y2, u2,v2, x3,y3, u3,v3];
% Tempi ai quali e' richiesta la soluzione.
ts = 0:0.004:5;
% Risoluzione del sistema con ode45 del Matlab.
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y] = ode45(@TremasseFun,ts,y0,options);
% Grafica.
figure(1);    
for it = 1:length(t)
    plot(y(it,1),y(it,2),'r.','markersize',16); hold on 
    plot(y(it,5),y(it,6),'w.','markersize',16)
    plot(y(it,9),y(it,10),'b.','markersize',16)
    plot(y(1:it,1),y(1:it,2),'r',y(1:it,5),y(1:it,6),'w',y(1:it,9),y(1:it,10),'b')
    title('Moto di tre masse'); axis([-1.5 2.5 -1.5 2.5]);   axis square
    hold off;
    drawnow
end
end




