clc; close all; clear;

filename = 'o2-ch4.txt';

% Leggi il file
txt = fileread(filename);

% Trova O/F
OF_cell = regexp(txt, 'O/F=\s*([-+]?\d*\.?\d+)', 'tokens');

% Trova T
T_cell = regexp(txt, 'T,\s*K\s+([-+]?\d*\.?\d+)', 'tokens');

% Converti le stringhe in numeri
OF = cellfun(@(x) str2double(x{1}), OF_cell);
T  = cellfun(@(x) str2double(x{1}), T_cell);

% Crea tabella
risultati = table(OF(:), T(:), ...
    'VariableNames', {'OF', 'T_K'});

disp(risultati);
plot(risultati, "OF", "T_K");
max(risultati.T_K)