close all; clc; clear all;

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
%  
% Codice di esempio che simula l'equazione di diffusione instazionaria con
% condizioni al contorno alla Dirichlet-Dirichlet (DD), Neumann-Dirichlet 
% (N-D) o periodiche (PER) mediante uno schema 'theta-method'.
% La equazione è:
% 
%                   phi_t = k*phi_xx + p
% 
% La codifica generale dello schema e' data da:
% 
%                  A*phi(n+1) = B*phi(n) + q
% 
% con A = I - beta*teta*D2 e B = I + beta*(1-teta)*D2.
% Le condizioni alla Dirichlet o Neumann sono codificate secondo la 
% notazione:
%       phi(0)    = PhiA
%       phi(L)    = PhiB
%    (dphi/dx)(0) = g0
% quelle periodiche sono:
%       phi(0) = phi(L)
% 
L     = 1;                    % Ampiezza del dominio (x in [0,L])
N     = 30;                   % Numero di punti interni al dominio
x     = linspace(0,L,N);      % Mesh spaziale
h     = x(2)-x(1);            % Passo spaziale (uniforme)
beta  = 1;                    % Assegnazione parametro beta = k*Dt/Dx^2
theta = 0.3;                  % Parametro numerico dello schema
k     = 1;                    % Coeff. di diffusione per il calcolo di Dt
Dt    = (beta*h^2)/k;         % Dt 'fisico'
p     = zeros(N,1);           % Produzione
T     = 1;                    % Tempo finale
Nt    = round(T/Dt);          % Numero di passi temporali

BC    = 'PER';                % Flag per la scelta delle cond. al contorno.

% Codifichiamo separatamente le diverse condizioni al contorno.
switch BC
    case 'D-D'                            % BCs alla Dirichlet-Dirichlet.
        Phi0   = -sin((3/2)*pi*x/(L))';   % Condizione iniziale.
% Costruzione dell'operatore 'base' di diffusione.
        D2     = -gallery('tridiag',N-2);
        I      = eye(N-2);
        A      = I - beta*theta*D2;
        B      = I + beta*(1-theta)*D2;
        PhiA   = 0;    PhiB   = 1;      % Valore delle cond. al contorno.
        q      = Dt*p(2:N-1);           % Termine noto nei punti interni.
        Phi    = Phi0(2:N-1);           % Cond. iniziale nei punti interni.        
        q(1)   = q(1)   + PhiA*beta;    % Modifica del termine noto per...
        q(end) = q(end) + PhiB*beta;    % ...tenere conto delle BCs.
        t      = 0;                     % Inizializzazione tempo.
        for it = 1:Nt
            t   = t + Dt;
            Phi = A\(B*Phi + q);            % Risoluzione
            plot(x,Phi0,'.-r',x,[PhiA;Phi;PhiB],'.-k')
            title(['Eq.ne Diffusione con BCs alla Dirichlet-Dirichlet.'])
            text(L/6,0.8,['it = ',num2str(it),'   t = ',num2str(t)]);
            axis([0 L -1 1]); drawnow
        end
    case 'N-D'                          % BCs alla Neumann-Dirichlet.
        Phi0   = cos(.5*pi*x/(L))';     % Condizione iniziale.
% Costruzione dell'operatore 'base' di diffusione.
        D2      = -gallery('tridiag',N-2);
% Modifica per la condizione alla Neumann in x=0
        D2(1,1) = D2(1,1) + 1;
        I       = eye(N-2);
        A       = I - beta*theta*D2;
        B       = I + beta*(1-theta)*D2;
        g0      = 1;    PhiB   = 0;      % Valore delle cond. al contorno.
        q       = Dt*p(2:N-1);           % Termine noto nei punti interni.
        Phi     = Phi0(2:N-1);           % Cond. iniziale nei punti interni.        
        q(1)    = q(1)   - g0*h*beta;    % Modifica del termine noto per...
        q(end)  = q(end) + PhiB*beta;    % ...tenere conto delle BCs.
        t       = 0;                     % Inizializzazione tempo.
        for it = 1:Nt
            t   = t + Dt;
            Phi = A\(B*Phi + q);         % Risoluzione
            plot(x,Phi0,'.-r',x,[Phi(1)-g0*h;Phi;PhiB],'.-k')
            title(['Eq.ne Diffusione con BCs alla Neumann-Dirichlet.'])
            text(2*L/3,0.8,['it = ',num2str(it),'   t = ',num2str(t)]);
            axis([0 L -1 1]); drawnow
        end        
    case 'PER'                           % BCs periodiche.
        Phi0  = sin(pi*x/(L))';          % Condizione iniziale.
% Costruzione dell'operatore 'base' di diffusione.
        v  = zeros(1,N-1);  v(1) = -2;  v(2) = 1;  v(end) = 1;
        D2 = gallery('circul',v);  
        I      = eye(N-1);
        A      = I - beta*theta*D2;
        B      = I + beta*(1-theta)*D2;
        p      = Dt*p(1:N-1);            % Termine noto nei punti interni.
        Phi    = Phi0(1:N-1);            % Cond. iniziale nei punti interni.        
        t      = 0;                      % Inizializzazione tempo.
        for it = 1:Nt
            t   = t + Dt;
            Phi = A\(B*Phi + p);         % Risoluzione
            plot(x,Phi0,'.-r',x,[Phi;Phi(1)],'.-k')
            title(['Eq.ne Diffusione con BCs Periodiche.'])
            text(L/3,0.2,['it = ',num2str(it),'   t = ',num2str(t)]);
            drawnow
        end
end
        
