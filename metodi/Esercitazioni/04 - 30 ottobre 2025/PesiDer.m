function W = PesiDer(xs,xc,p)

% Esercitazione di 'Metodi Numerici in Ingegneria Aerospaziale'
% Prof. G. Coppola
% 
% Function per il calcolo dei pesi delle formule di derivazione mediante il
% metodo dei coefficienti indeterminati.
%
% Parametri di ingresso:
% xs   -->    Vettore di coordinate dei nodi del mesh.
% xc   -->    Punto di collocazione delle derivate.
% p    -->    Ordine di derivazione della formula richiesta.
% 
% Parametri di uscita:
% W    -->    Riga dei pesi delle derivate o, se nargin = 2, matrice
%             contenente sulla j-ma riga i pesi della derivata j-1.

N   = length(xs);          % N è il numero di nodi dello stencil.
xs  = xs(:).';             % Rendiamo il vettore x una riga. 
csi = xs - xc;             % Definiamo la variabile ausiliaria x-xc.
M   = zeros(N);            % Allocazione della matrice M.
% Composizione della matrice.
for j = 0:N-1
    M (j+1,:) = (csi.^j)/factorial(j);
end
I = eye(N);               % Matrice identica.
% Risoluzione del sistema per i pesi di tutte le derivate (se è fornito p, 
% il termine noto per i pesi della derivata di ordine p sarà la (p+1)-ma 
% colonna della matrice identica).
W = M\I;                  % Risoluzione del sistema.
W = W.';                  % Mettiamo i pesi sulle righe.
% Se richiesto, estraiamo solo i pesi della p-ma derivata (NB la presente
% versione non è 'ottimizzata', poichè calcola in ogni caso i pesi di tutte
% le derivate, anche se è richiesta solo la derivata di ordine p).
if nargin==3
    W = W(p+1,:);
end
