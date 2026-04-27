function k = kappa(x)

global L

% Legge di variazione del coefficiente di diffusione k(x).

k = 1 + 2*(x/L).^2;