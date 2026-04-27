function dPhidt = RHS(~,Phi)

global Op 

Phi = Phi(:);
dPhidt = Op*Phi;

end