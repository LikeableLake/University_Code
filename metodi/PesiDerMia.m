function W = PesiDer(xs,xc,d)

%routine per il calcolo dei pesi nelle derivate formule backward, forward e
%centrale
%formula di derivazione f^(n) = sum((c_k)*(f_k)) k=1:n
%dagli sviluppi di taylor f_k=sum(c_k) k=1:inf * sum(f^(j) (x_c)* ((x_k-x_c)^j)/j! j=0:inf
%=sum(f^(j) (x_c))j=0:inf * sum(c_k ((x_k -x_c)^j)/j!)k=1:n

%xs=elenco punti
%xc=punto in cui derivare
%d=grado derivata

%devo scrivere la matrice dei c_k

%1      1       1       1...
%x1-xc  x2-xc   x3-xc   ...
%...
%(x1-xc)^j /j!...

%num2str (123) restituisce 123 come stringa

xs  =   xs(:)';   %trasforma xs in sicuramente un vettore riga

csi =   xs-xc; %a ogni componente del vettore xs(riga) va sottratto lo scalare xc
n   =   length(xs);
M   =   zeros(n);

for j=1:n
    M(j,:)=(csi.^(j-1))/ factorial(j-1);    %la riga della matrice colonna è csi elevata a un j, ovvero elevo tuttle le sue componenti alla j
end


%inv di A è pesante per risolvere i sistemi e poco precisa, matlab risolve con:
%inv(A)*b= A\b
%b*inv(A)=b/A
%sono algoritmi interni migliori
%inv(A)=A\I

%uso nargin per capire quanti dati mi sono stati dati
    I=eye(n);

if nargin==3
    b=I(:,d+1); %d-esima colonna della matrice identica

    W=M\b;      %d-esima colonna della matrice M ovvero derivata d-esima
else
            %oppure per trovaare tutti i pesi
    W=M\I;  %fornisce tutti i pesi come colonne di W
end
%in questo modo io posso chiedere, usando solo xc e xs, la matrice
%completa, oppure aggiungendo il grado di derivazione, la colonna precisa
W=W';

end