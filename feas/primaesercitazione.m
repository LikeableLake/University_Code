clc; clear; close all;

c=3*10^8;       %velocita della luce 
f=100e6;        %frequenza

lambda=c/f;

d=dipole('width', lambda/100, 'length', lambda/10);
%costruisce un dipolo elementare e ne inserisce tutti i dati in una classe


show(d);

%il dipolo non e di forma cilindrica ma schiacciata perche lo si vuole
%mettere di default in un'antenna
%il comportamento, nonostante non e perfettamente cilindrico sara molto
%prossimo poiche la sezione e molto piccola rispetto alla lunghezza d'onda

%feed = entrata
%pec = perfect electric conductor

pattern(d, f,);

%diagramma di radiazione del dipolo(toroidale) senza espansioni polari
%quindi in realta e un po' diverso da quello che abbiamo studiato, irradia
%alla stessa maniera (va dimostrato)

%il diagramma è in scala logaritmica (dB)
