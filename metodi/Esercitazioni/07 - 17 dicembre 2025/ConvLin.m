close all; clc; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che simula la equazione di convezione lineare con 
% condizioni al contorno periodiche.
% La equazione e' discretizzata con vari schemi classici che possono essere
% selezionati tramite la variabile 'Schema'.
% 
% 'CEN'        Schema FTCS
% 'UPW'        Schema upwind (per a>0)
% 'DWN'        Schema downwind (per a>0)
% 'LF'         Schema di Lax-Friedrichs
% 'LW'         Schema di Lax-Wendroff 
% 'BW'         Schema di Beam-Warming 
% 'LEAP'       Schema Leapfrog
%

L  = 1;                     % Ampiezza del dominio (x in [0,L]).
N  = 100;                   % Numero di punti interni al dominio.
x  = linspace(0,L,N);       % Mesh spaziale.
h  = x(2)-x(1);             % Passo spaziale (uniforme).
a  = 1;                     % Velocità di convezione.
C  = 0.5;                   % Numero di Courant C=a*Dt/h
Dt = C*h/abs(a);            % Step temporale.
T  = L/abs(a);              % Tempo finale (1 giro)
% Ricalcoliamo il Dt in modo tale che in Nt*Dt compio un 'giro'.
Nt = round(T/Dt);   Dt = T/Nt;   C = a*Dt/h;

IC = 3;                     % Selezione della condizione iniziale
switch IC
    case 1                  % Gaussiana
        s = 100/L^2;   phi0= exp(-s*(x-L/2).^2);
    case 2                  % Onda sinusoidale
        phi0= sin(2*pi*x/L);
    case 3                  % Onda quadra
        phi0 = zeros(N,1);
        for i= 1:N
            if (x(i)>L/3) && (x(i)<2*L/3)
                phi0(i) = 1;
            end
        end
end
% Set dello schema numerico 
SchemaAll ={'CEN','UPW','LF','LW','BW','LEAP'};
for caso = 1:6
    Schema = SchemaAll{caso};
% Set del valore al tempo t=0.
    Phi = phi0(1:N-1);    Phi = Phi(:);
% Costruzione della matrice di transizione
    switch Schema
       case 'CEN'               % Schema FTCS
            v = zeros(1,N-1); v(1) = 1; v(2) = -C/2; v(end) = C/2; 
            A = gallery('circul',v);
            it = 0;              tit = 'FTCS';                
        case 'UPW'              % Schema upwind
            v = zeros(1,N-1);    v(1) = 1-C;            v(end) = C;
            A = gallery('circul',v);
            it = 0;              tit = 'Upwind (per a>0)';
        case 'DWN'              % Schema downwind
            v = zeros(1,N-1);    v(1) = 1+C;            v(2) = -C;
            A = gallery('circul',v);
            it = 0;              tit = 'Downwind (per a>0)';
        case 'LF'               % Schema di Lax-Friedrichs
            v = zeros (1,N-1);   v(2)   = (1-C)/2;      v(end) = (1+C)/2;
            A = gallery('circul',v);
            it = 0;              tit = 'Lax-Friedrichs';
        case 'LW'               % Schema di Lax-Wendroff
            v = zeros (1,N-1);   v(1) = (1-C^2);        v(end) = C*(C+1)/2;
            v(2) = C*(C-1)/2;
            A = gallery('circul',v);
            it = 0;              tit = 'Lax-Wendroff';
        case 'BW'               % Schema di Beam-Warming
            v = zeros (1,N-1);   v(1) = (2-3*C+C^2)/2;  v(end) = -C*(C-2);
            v(end-1) = C*(C-1)/2 ;
            A = gallery('circul',v);
            it = 0;              tit = 'Beam-Warming';
        case 'LEAP'              % Schema Leapfrog
            % Primo step: Lax-Wendroff
            v = zeros (1,N-1);   v(1) = (1-C^2);        v(end) = C*(C+1)/2;
            v(2) = C*(C-1)/2;
            A = gallery('circul',v);
            PhiOld = Phi;        Phi = A*PhiOld;
            v = zeros (1,N-1);   v(2) = -C;             v(end) = C;
            A = gallery('circul',v);
            it = 1;              tit = 'Leapfrog';
    end
% Ciclo instazionario
    while it < Nt
        it = it +1;
        if strcmp(Schema,'LEAP')
            PhiNew = A*Phi + PhiOld;
            PhiOld = Phi; Phi = PhiNew;
        else
            Phi = A*Phi;
        end
% Grafica
        if mod(it,5)==0
            subplot(2,3,caso);
            plot(x,phi0,'r.-',x,[Phi;Phi(1)],'k.-');
            axis([0 L -0.5 1.5]);     title(tit);
            drawnow
        end
    end 
end
