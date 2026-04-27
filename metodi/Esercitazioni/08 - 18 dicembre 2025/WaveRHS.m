function dPhidt = WaveRHS(~,Phi)

% Function utilizzata dall'integratore ODE per la equazione delle
% onde con varie BCs 

global A

Phi    = Phi(:);
dPhidt = A*Phi;

end