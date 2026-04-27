function WaveMain(caso)

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Codice di esempio che simula la equazione delle onde con condizioni
% al contorno omogenee (D-D, D-N, ecc... a seconda del valore dell'unico 
% input 'caso') mediante una tecnica di semi-discretizzazione e l'utilizzo 
% di integratori ODE del Matlab.
% 
% Il codice riduce il sistema di N equazioni differenziali ordinarie
% del secondo ordine, cui si perviene tramite una discretizzazione
% spaziale dell'operatore di derivata II, in un sistema di 2*N
% equazioni differenziali del I ordine, mediante la definizione di
% nuove variabili pari alle derivate delle variabili originarie.
% Con tale trasformazione il sistema di ODE originario:
%
%                        y_tt = cq*D2*y
% 
% viene ricondotto al sistema di ODE del I ordine:
% 
%                        _ _   _                   _   _ _
%                       |   | |                     | |   |
%                    d  |y_1| |      O         I    | |y_1|
%                   ----|   |=|                     |.|   |
%                    dt |y_2| | ((c/dx)^2)*D2  O    | |y_2|
%                       |_ _| |_                   _| |_ _|
%

global A % Dichiaraz. di variabili globali da passare 
         % alla function per la integrazione della ODE.  

% Costanti e parametri fisici                              
L = 1;                        % Lunghezza del dominio
c = 1;   cq = c*c;            % Velocita'
T = 2*L/c;                    % Tempo finale

% Discretizzazione
N  = 100;                     % Numero di nodi spaziali
x  = linspace(0,L,N);         % Mesh spaziale
h  = x(2) - x(1);  hq = h*h;  % Passo spaziale
x  = x(2:end-1);              % Eliminazione variabili di bordo
Nt = 100;                     % Numero di nodi temporali
t  = linspace(0,T,Nt);        % Mesh temporale

% Condizioni iniziali
% La PDE e' del II ordine nel tempo: sono richiste condizioni iniziali sia 
% sulla posizione che sulla velocita'.

switch caso % Vari casi. In dipendenza da questi va modificato anche
            % l'operatore differenziale D2 (e quindi A)
            
    case 1  % BC D-D
            % Impulso gaussiano con velocita' iniziale nulla, cosi' che si
            % propaghi in due versi, i.e. lungo entrambe le caratteristiche
        y0  = .1*exp(-100*(x-L/2).^2);
        yp0 = 0*y0;

    case 2  % BC D-D
            % Impulso gaussiano con velocita' iniziale tale da farlo propagare
            % in un solo verso, i.e. lungo una sola caratteristica
        y0  = .1*exp(-100*(x-L/2).^2);      % Valore iniziale
        yp0 = -c*(-100*2*(x-L/2)).*y0;      % Pendenza iniziale
        
    case 3  % BC D-D
            % Onda stazionaria
        y0  = .1*sin(pi*x/L); 
        yp0 = 0*y0;
        
    case 4  % BC D-N
            % Impulso gaussiano con velocita' iniziale nulla
        y0  = .1*exp(-100*(x-L/2).^2);
        yp0 = 0*y0;
        
    case 5  % per BC N-N
            % Impulso gaussiano con velocita' iniziale nulla
        y0  = .1*exp(-100*(x-L/2).^2);
        yp0 = 0*y0;
end

% Operatori differenziali
% Costruzione dell'operatore a blocchi
O  = zeros(N-2);
I  = speye(N-2);
D2 = gallery('tridiag',N-2,1,-2,1); % BC alla Dirichlet omogenee

% Implementazione delle BC
% Per le condizioni iniziali scelte con lo switch precedente sono
% previsto modifiche dell'operatore differenziale in maniera da
% implementare condizioni al contorno coerenti
switch caso
    case {1,2,3} % BC D-D
        % nessuna modifica
        
    case 4 % BC D-N
        D2(end,end) = D2(end,end) + 1;
        
    case 5 % BC N-N
        D2(1,1)     = D2(1,1)     + 1;
        D2(end,end) = D2(end,end) + 1;
end
D2 = D2/hq;
% Assembliamo l'operatore relativo al sistema di ODE del primo ordine
A  = [  O    I 
      cq*D2  O ];

% Integrazione temporale
% Assegnazione tolleranze per l'integratore ODE.
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
% Integraziome della ODE.
[t,y]   = ode45(@WaveRHS,t,[y0 yp0],options);

% Grafica instazionaria
massimo = max(max(abs(y(:,1:N-2))));
canvas = [0 L -1.2*massimo 1.2*massimo];
for i = 1:Nt
    figure(1)
    plot(x,y0,'.-r',x,y(i,1:N-2),'.-k')
    title('Wave equation')
    axis(canvas), grid on
    drawnow
end

% Grafica finale.
figure(2)
colormap('gray')
waterfall(x,t,y(:,1:N-2),ones(size(y(:,1:N-2))));
view(10,70)

end
