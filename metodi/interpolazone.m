% interpolazione polinomiale approccio van der monde
% P(x)= a0 + a1*x + a2*x^2...
% P(x_i)=phi_i
% la matrice del sistema è:
% 1 x_0 x_0^2 x_0^3     a_0     phi_0
% 1 x_1 x_1^2 x_1^3  *  a_1  =  phi_1
% 1 x_2 x_2^2 x_2^3     a_2     phi_2


clear; close all; clc

N=10;        %grado del polinomio
g=N-1;      
L=1;        %valore finale dello stencil
x= linspace(0,L,N);     %stencil sotto forma di vettore
x=x(:);

V=zeros(N);

for i= 1:N
    V(:,i)= x.^(i-1);       %matrice di van der monde
end

f= x.*cos(pi*x/L);      %valori della funzione da interpolare nei corrispondenti x

a=V\f;      %coefficienti del polinomio

%creo un mesh molto più fitto in modo da rappresentare in grafico il
%polinomio interpolante

Ng=100;
xg=linspace (-1,L+1,Ng);
fg = xg.*cos(pi*xg/L);

P=0;
for ig=0:g
    P= P+ a(ig+1)*xg.^ig;       %costruzione dei valori del polinomio
end

plot (xg, fg, '-r'); axis ([-1 L+1 -2 1]); hold on;

plot (x, f, 'ow', xg, P, '-b');



%studia comando flip, polyfit(che sarebbe tutto il bordello scritto finora
%in un singolo comando)+polyval (che invece restituisce valori)

%fg=polyval(polyfit(x, f, g), xg)



%% interpolazione con lagrange
close all; clear; clc;

N=5;     %grqado del polinomio
g=N-1;
L=1;
x=linspace(0, L, N);
I=eye(N);
Ng=100;
xg = linspace(-L/2, 1.5*L, Ng);

for i=1:N
    f= I(:,i);
    fg= polyval(polyfit( x, f, g), xg);
    plot (xg, fg, '-r', x, f,'o'); axis ([-L/2 1.5*L , -1 2]); hold on;
    pause
end


%% interpolazione della curva di Runge
close all; clear; clc;

N=15;
Ng=200;
g=N-1;
L=4;
x=linspace(-L, L, N);
f= 1./(1+x.^2);
plot(x,f, 'o'); axis([-L L, -1 2]); hold on;
xg=linspace(-L, L, Ng);
fg= 1./(1+xg.^2);
plot(xg,fg); hold on;
Fl= polyval(polyfit(x,f,g),xg);
plot(xg,Fl,'r-');

%comando subpot  mette più grafici in una matrice grafica

%% mesh di cerbiscev

teta=linspace(pi, 0, N);
xc=L*cos(teta);

% figure(2);
% plot(xc, 0*xc, 'or', x, 1+0*x, 'ow'); axis([-L L, -1 2]);

fc= 1./(1+xc.^2);
Fc= polyval(polyfit(xc,fc,g),xg);
plot( xg, Fc);


figure(2);
subplot(3,1,1); plot(x,f, 'o'); hold on; axis([-L L, 0 1.2]); plot(xg,fg); hold on; plot(xg,Fl,'r-');
subplot(3,1,2); plot(x,f, 'o'); hold on; axis([-L L, 0 1.2]); plot(xg,fg); hold on; plot( xg, Fc);

%% spline

ps=spline(x, f, xg); %spline che passa per i valori f sul mesh x stimata su xg punti
subplot(3,1,3); plot(x,f, 'o'); hold on; axis([-L L, 0 1.2]); plot(xg,fg); hold on; plot( xg,ps);

%la spline interpola megliuo di tutte perchè è un insieme di cubiche quindi
%non oscilla come invece i polinomi

%se non metto il terzo valore spline mi fornisce dei valori da usare in
%ppval (una piecewise polynomial structure), un array strutturato
% 
% pp=spline(x,f);
% pp
%PP=ppval(pp, xg)
%la curva ha come sistema di riferimento lo 0, quindi per metterla nel mio
%plot al posto di usare ppval posso estrarre ciò che mi interessa dalla
%struct   pp.coefs restituisce la matrice con i coefficienti, 4
%coefficienti per n cubiche che servono



%%