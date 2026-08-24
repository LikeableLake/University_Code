clc; clear; close all

filename = 'o2-ch4-50atm.txt';

% Leggi il file
txt = fileread(filename);

% Trova O/F
OF_cell = regexp(txt, 'O/F=\s*([-+]?\d*\.?\d+)', 'tokens');

% Trova T
T_cell = regexp(txt, 'T,\s*K\s+([-+]?\d*\.?\d+)', 'tokens');
CO2_cell = regexp(txt, '*CO2\s+([-+]?\d*\.?\d+)', 'tokens');

% Converti le stringhe in numeri
OF = cellfun(@(x) str2double(x{1}), OF_cell);
T  = cellfun(@(x) str2double(x{1}), T_cell);
CO2 = cellfun(@(x) str2double(x{1}), CO2_cell);
% Crea tabella
risultati = table(OF(:), T(:), CO2(:), 'VariableNames', {'OF', 'T_K', 'CO2'});

disp(risultati);
yyaxis left; plot (risultati, "OF", "T_K"); hold on; 
yyaxis right;
plot (risultati, "OF", "CO2")


max(risultati.T_K)