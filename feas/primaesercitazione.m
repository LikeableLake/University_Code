clc; clear; close all;

c=3*10^8;       %velocita della luce 
f=100e6;        %frequenza

lambda=c/f;

d=dipole('width', lambda/100, 'length', lambda/10);
%costruisce un dipolo elementare e ne inserisce tutti i dati in una classe


show(d);

%il dipolo non è di forma cilindrica ma schiacciata perche lo si vuole
%mettere di default in un'antenna
%il comportamento, nonostante non è perfettamente cilindrico sara molto
%prossimo poiche la sezione e molto piccola rispetto alla lunghezza d'onda

%feed = entrata
%pec = perfect electric conductor

figure(1);
pattern(d, f);

%diagramma di radiazione del dipolo(toroidale) senza espansioni polari
%quindi in realtà è un po' diverso da quello che abbiamo studiato, irradia
%alla stessa maniera (va dimostrato)

%il diagramma è in scala logaritmica (dB)

%diagramma di direttività, ossia direzione per direzione esprime il campo,
%in questo caso in dB quindi scalato di una certa quantita

figure(2);
pattern(d, f, 'type', 'efield');

%con questa formula il campo è in scala lineare e in volt su metro, in
%questo modo la forma è più simile al toro che ci aspettavamo

%come mai si usano i dB? solitamente contraggono le scale è utile in
%circostanze in cui sono presenti contemporaneamente lobi molto grandi e
%lobi molto piccoli, come nelle antenne molto direzionali, che si unsano
%per evitare di disturbare antenne che non devono ricevere tale segnale

%il dipolo non ha lobi laterali diversi dai suoi due lobi principali

%di solito i lobi laterali di un'antenna direzionale saranno 20 dB sotto il
%lobo principale in intensità

%di solito non si dà mai il solido di rotazione ma due tagli su piani
%tipicamente ortogonali, nel caso del dipolo ad esempio un taglio che
%contenga l'asse z ossia il disegno dei lobi e un altro sul piano xy che
%conterrà un cerchio (senza il suo centro?)

%il comando patternelevation fa vedere il pattern in funzione
%dell'elevazione, ossia l'opposto di theta

%il campo lontano ha una banda finita e si può applicare il teorema del
%campionamento per definire con quale passo angolare campionare il campo

%in particolare quanto varia angolarmente il campo lontano è legato alla
%dimensione dell'antenna e in particolare alla trasformata di fourier
%dell'antenna


%__________________________________________________________________________

% dir=patternElevation(d, f, 0, 'Elevation', (-90:1:90)); 
% figure(3)
% 
% polarplot(dir)

%_______________________questo pezzo non funziona__________________________

figure(3)
subplot(1,2,1); dir=pattern(d,f,0,(-180:1:180), 'Type','efield');
polarplot(dir); title('direttività piano xz')

%al posto di patternElevation posso usare pattern e fissare l'azimuth al
%variare dell'elevation

subplot(1,2,2); az=pattern(d,f,(-180:1:180), 90, 'type','efield');
polarplot(az*1e9); title('direttività piano xy (moltiplicata per 10^9)')

%qui fisso l'elevation (a 90 poichè serve l'angolo theta) e faccio variare 
%%l'azimuth